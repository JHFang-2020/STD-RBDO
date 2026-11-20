function [Snew,info] = c36_acq_sur(model,bank,evs,dw,B,mode,opts)
%C36_ACQ_SUR  Batch stepwise uncertainty reduction over spatiotemporal events.
%   evs  : cell of event structs (c36_events), one per design in the design
%          neighbourhood; dw : their weights (sum 1)
%   mode : 'closed' - ST-MSUR+: trajectory-event SUR, closed form (Phi2) under
%                     the Clark margin, greedy batch with exact Schur updates
%          'mc'     - ST-MSUR (v35 as written): capped expected-U screening of
%                     the event points, shortlist re-ranked by Monte-Carlo
%                     trajectory-event SUR (32 fantasies x 128 conditional draws)
%          'point'  - pointwise SUR ablation: same closed form, but every
%                     population point (n,t) is its own excursion event
%   A snapshot at state c observes all spatial components z_l(c) exactly,
%   i.e. a complete retained observation block.
if nargin<7, opts = struct(); end
nExplore = getf(opts,'nExplore',200);
seed     = getf(opts,'seed',1);
Ucap     = getf(opts,'Ucap',6);
nShort   = getf(opts,'shortlist',8);
nFan     = getf(opts,'nFantasy',32);
nCond    = getf(opts,'nConditional',128);
c = c36_config(); N = bank.N; r = model.r;
rs = RandStream('mt19937ar','Seed',seed);
% ---------------- assemble events (E) and their points (X)
Ex = {}; Ej = {}; Ee = {}; Ew = {}; Pm = {}; Pv = {};
Em = []; Evv = []; EW = []; Ek = {}; nEv = 0; Cev = {};
for j = 1:numel(evs)
    ev = evs{j}; if isempty(ev), continue; end
    if strcmp(mode,'point')
        lin = find(~isnan(ev.sig));
        pp = 0.5*erfc(ev.G(lin)./(sqrt(2)*max(ev.sig(lin),1e-12)));
        keep = pp.*(1-pp) > 1e-7; lin = lin(keep);
        [~,o] = sort(pp(keep).*(1-pp(keep)),'descend'); lin = lin(o(1:min(end,3000)));
        for q = 1:numel(lin)
            nEv = nEv+1; Ex{nEv} = lin(q); Ej{nEv} = evs{j}.Pj(lin(q)); Ew{nEv} = 1;
            Em(nEv,1) = ev.G(lin(q)); Evv(nEv,1) = ev.sig(lin(q))^2; EW(nEv,1) = dw(j)/bank.nt; %#ok<AGROW>
        end
    else
        for e = 1:numel(ev.n)
            k = ev.npts(e); nEv = nEv+1;
            Ex{nEv} = ev.pts(e,1:k)'; Ej{nEv} = evs{j}.Pj(ev.pts(e,1:k)'); Ew{nEv} = ev.w{e};
            Em(nEv,1) = ev.m(e); Evv(nEv,1) = ev.v(e); EW(nEv,1) = dw(j); %#ok<AGROW>
            Pm{nEv} = ev.mu{e}; Cev{nEv} = ev.C{e};
        end
    end
end
info = struct('mode',mode,'nEvents',nEv);
if nEv == 0
    Snew = maximin(bank,model.S,B,rs); info.fallback = 'maximin (no ambiguous event)'; return;
end
Evv = max(Evv,1e-300); aE = -Em./sqrt(Evv);
% points: unique (state, node) pairs
allLin = vertcat(Ex{:}); allJ = vertcat(Ej{:});
[uState,~,mapS] = unique(allLin);
nX = numel(allLin);
evOf = repelem((1:nEv)',cellfun(@numel,Ex));
wOf  = vertcat(Ew{:});
ptEnd = cumsum(cellfun(@numel,Ex))'; ptStart = [1; ptEnd(1:end-1)+1];
% ---------------- candidates
expl = rs.randperm(N*bank.nt,min(nExplore,N*bank.nt))';
candLin = unique([uState; expl]);
Sc = bank.S(candLin,:); nC = numel(candLin);
[Wc,zvC] = c36_cw(model,Sc);
Sx = bank.S(uState,:);
[Wx,~] = c36_cw(model,Sx);
Xe = cell(1,r); Xp = cell(1,r);
needPts = strcmp(mode,'mc');
for l = 1:r
    cm = model.comp{l};
    Cc = c36_ckern(cm,Sx,Sc) - Wx{l}'*Wc{l};                            % nU-by-nC
    rowsPt = -model.Phi(allJ,l).*Cc(mapS,:);                            % nX-by-nC (cov of g with z_l)
    A = sparse(evOf,(1:nX)',wOf,nEv,nX);
    Xe{l} = full(A*rowsPt);                                             % nEv-by-nC
    if needPts, Xp{l} = rowsPt; end
end
Vc = zvC;                                   % r-by-nC conditional variances
s2E = zeros(nEv,1);
selected = zeros(B,1); scores = zeros(B,1);
Kcs = cell(1,r); for l = 1:r, Kcs{l} = zeros(nC,0); end
if needPts
    ptMu = zeros(nX,1); ptV = zeros(nX,1); pos = 0;
    for e = 1:nEv
        k = numel(Ex{e}); ptMu(pos+(1:k)) = Pm{e}; ptV(pos+(1:k)) = diag(Cev{e}); pos = pos+k;
    end
    ptV = max(ptV,1e-300); s2P = zeros(nX,1);
    ptW = EW(evOf)./cellfun(@numel,Ex(evOf))';
    Bphi = cell(nEv,1); for e = 1:nEv, Bphi{e} = zeros(numel(Ex{e}),0); end
    rsMC = RandStream('mt19937ar','Seed',seed+77);
    Zf = rsMC.randn(nFan, r*B);                         % common fantasy normals
    Zc = cell(nEv,1); Zb = cell(nEv,1);
    for e = 1:nEv
        Zc{e} = rsMC.randn(numel(Ex{e}),nCond*nFan);    % common conditional normals
        Zb{e} = rsMC.randn(numel(Ex{e}),nCond);
    end
    Hcur = 0;
    for e = 1:nEv
        Rt = psdroot(Cev{e});
        pe = mean(any(Pm{e} + Rt*Zb{e} <= 0,1));
        Hcur = Hcur + EW(e)*nCond/(nCond-1)*pe*(1-pe);
    end
end
if strcmp(mode,'tmp')
    [~,ptCand] = ismember(allLin,candLin);
end
for b = 1:B
    valid = true(nC,1);
    if b>1, valid(selected(1:b-1)) = false; end
    if strcmp(mode,'tmp')
        % TMP-greedy (trajectory misclassification probability, cf. AKFDR,
        % Song et al. 2026): most ambiguous trajectory event under the
        % Kriging-believer variance, sampled at its most uncertain critical point
        vRem = max(Evv - s2E,1e-300);
        pK = 0.5*erfc(Em./sqrt(2*vRem));
        tmpv = EW.*min(pK,1-pK);
        [~,eo] = sort(tmpv,'descend'); cs = [];
        for ee = eo'
            pr = ptStart(ee):ptEnd(ee);
            sdk = sqrt(max(diag(Cev{ee}),0)).*abs(Ew{ee});
            [~,ko] = sort(sdk,'descend');
            for kk = ko'
                ci = ptCand(pr(kk));
                if ci>0 && valid(ci), cs = ci; break; end
            end
            if ~isempty(cs), scores(b) = tmpv(ee); break; end
        end
        if isempty(cs), cs = find(valid,1); end
    elseif needPts
        % --- screening: capped expected future U of the event points
        ds2 = zeros(nX,nC);
        for l = 1:r
            ok = Vc(l,:) > 1e-12*model.comp{l}.sf2;
            ds2(:,ok) = ds2(:,ok) + Xp{l}(:,ok).^2./Vc(l,ok);
        end
        expl2 = min(s2P + ds2, ptV*(1-1e-12));
        EU = futureU(ptMu,ptV,expl2);
        U0 = abs(ptMu)./sqrt(ptV);
        gainS = ptW'*(min(EU,Ucap) - min(U0,Ucap));
        gainS(~valid) = -Inf;
        [~,o] = sort(gainS,'descend'); short = o(1:min(nShort,sum(valid)));
        % --- re-rank by Monte-Carlo trajectory-event SUR (CRN across candidates)
        gE = zeros(numel(short),1);
        for q = 1:numel(short)
            cq = short(q); fr = zeros(nFan,1);
            for e = 1:nEv
                pr = ptStart(e):ptEnd(e); k = numel(pr);
                Bn = zeros(k,r);
                for l = 1:r
                    if Vc(l,cq) > 1e-12*model.comp{l}.sf2
                        Bn(:,l) = Xp{l}(pr,cq)/sqrt(Vc(l,cq));
                    end
                end
                Bfull = [Bphi{e}, Bn];
                Cn = Cev{e} - Bfull*Bfull';
                Rt = psdroot((Cn+Cn')/2);
                shift = Bfull*Zf(:,1:size(Bfull,2))';                 % k-by-nFan
                F = Rt*Zc{e} + repelem(shift,1,nCond) + Pm{e};        % k-by-(nCond*nFan)
                pe = mean(reshape(any(F <= 0,1),nCond,nFan),1)';
                fr = fr + EW(e)*nCond/(nCond-1)*pe.*(1-pe);
            end
            gE(q) = Hcur - mean(fr);
        end
        [scores(b),qb] = max(gE); cs = short(qb);
    else
        % --- closed-form event SUR
        ds2 = zeros(nEv,nC);
        for l = 1:r
            ok = Vc(l,:) > 1e-12*model.comp{l}.sf2;
            ds2(:,ok) = ds2(:,ok) + Xe{l}(:,ok).^2./Vc(l,ok);
        end
        rho = min((s2E + ds2)./Evv, 1-1e-12);
        base = c36_phi2gain(aE,s2E./Evv);
        G = EW'*(c36_phi2gain(aE,rho) - base);
        G(~valid) = -Inf;
        [scores(b),cs] = max(G);
    end
    selected(b) = cs;
    % --- Schur update: condition (variances only) on the selected block
    for l = 1:r
        vstar = Vc(l,cs);
        if vstar <= 1e-12*model.comp{l}.sf2, continue; end
        cm = model.comp{l};
        kc = c36_ckern(cm,Sc,Sc(cs,:)) - Wc{l}'*Wc{l}(:,cs);
        if ~isempty(Kcs{l}), kc = kc - Kcs{l}*Kcs{l}(cs,:)'; end
        phiC = kc/sqrt(vstar);
        phiE = Xe{l}(:,cs)/sqrt(vstar);
        Vc(l,:) = max(Vc(l,:) - (phiC').^2, 0);
        Xe{l} = Xe{l} - phiE*phiC';
        s2E = s2E + phiE.^2;
        Kcs{l} = [Kcs{l}, phiC];
        if needPts
            phiP = Xp{l}(:,cs)/sqrt(vstar);
            Xp{l} = Xp{l} - phiP*phiC';
            s2P = s2P + phiP.^2;
            for e = 1:nEv, Bphi{e} = [Bphi{e}, phiP(ptStart(e):ptEnd(e))]; end
        end
    end
end
Snew = Sc(selected,:);
info.scores = scores; info.selectedLin = candLin(selected); info.nCandidates = nC;
info.explainedEventFraction = sum(EW.*min(s2E./Evv,1))/sum(EW);
end

function EU = futureU(mu,v,s2)
% E[|mu+|]/sqrt(v-s2), mu+ ~ N(mu,s2) (folded normal), cf. manuscript Eq. (22)
sv = sqrt(max(s2,0)); a = abs(mu);
z = a./max(sv,1e-300);
fold = sv.*sqrt(2/pi).*exp(-0.5*z.^2) + a.*erf(z/sqrt(2));
EU = fold./sqrt(max(v - s2,1e-300));
U0m = repmat(a./sqrt(v),1,size(s2,2));          % no information: future U = current U
EU(s2<=0) = U0m(s2<=0);
end

function R = psdroot(C)
[V,D] = eig((C+C')/2,'vector'); D = max(D,0);
R = V.*sqrt(D)';
end

function S = maximin(bank,Str,B,rs)
cand = rs.randperm(bank.N*bank.nt,min(4000,bank.N*bank.nt))';
X = bank.S(cand,:)./[6 6 4]; T = Str./[6 6 4];
d = min(pdist2(X,T),[],2); S = zeros(B,3);
for b = 1:B
    [~,k] = max(d); S(b,:) = bank.S(cand(k),:);
    d = min(d,pdist2(X,X(k,:)));
end
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
