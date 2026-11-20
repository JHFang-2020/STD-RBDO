function [Snew,info] = c36_acq_point(model,bank,ev,B,method,opts)
%C36_ACQ_POINT  Literature batch learning functions on the shared population.
%   All methods use the same surrogate, population, candidate points and the
%   same snapshot observation block; only the learning function / batch rule
%   differs.  Candidate points are the population points (n,t) at the spatial
%   minimiser p*, at the current design, with the common mode sampled.
%   'ukb'   : U = |mu|/sigma (Echard et al. 2011, AK-MCS) + Kriging Believer
%             batch (Ginsbourger et al. 2010: variance conditioned on the
%             selected blocks, mean frozen)
%   'silk'  : trajectory-level U of SILK (Hu & Mahadevan 2016, JMD 138:061406):
%             a trajectory is
%             resolved as failed once any point has mu<0 and U>=2; otherwise
%             its value is min_t U; batch by Kriging Believer
%   'efflp' : EFF (Bichon et al. 2008, eps = 2 sigma) with local penalisation
%             (Gonzalez et al. 2016; Wang, Monaghan & Zhang 2025):
%             phi_j(x) = P(|G_j| < L ||x-x_j||), L = max gradient norm
%   'pbalc' : learning function 3 of PBALC (Dang et al. 2024, Eq. 40):
%             [Phi(b-U) - Phi(-b-U)] f(x), b = 1, f = input density;
%             batch by Kriging Believer
%   'epak'  : EPAK (Hu, Dang, Wang, Beer & Wang 2025, RESS 262:111194):
%             trajectory misclassification p_w = Phi(-min_t U), batch by the
%             influence function IF = prod(1-R) times p_w (training + in-batch
%             points), snapshot at the trajectory's min-U time
%   'random': maximin space filling over population states (control)
if nargin<6, opts = struct(); end
seed = getf(opts,'seed',1); Mtop = getf(opts,'Mtop',3000);
rs = RandStream('mt19937ar','Seed',seed);
c = c36_config(); N = bank.N; nt = bank.nt; r = model.r;
info = struct('method',method);
if strcmp(method,'random')
    Snew = maximin(bank,model.S,B,rs); return;
end
if strcmp(method,'slabalq')
    % SL-ABALQ (Dang, Li, Valdebenito & Faes 2025, MSSP 241:113473):
    % time step (their Eq. 32): j* = argmax_j (1/N) sum_n s_nj^2 Phi(-r_n^min) Phi(r_n^min),
    %   with the trajectory-level minimum ratio r_n^min = min_j mu_nj/s_nj;
    % sample step (Eqs. 30/34, density factor dropped for MC candidates):
    %   n* = argmax_n s^2 Phi(-r) Phi(r) at time j*; batch completed by Kriging Believer
    sig = ev.sig; sig(isnan(sig)) = 0;
    linAll = find(sig>0);
    [Wx,Vx] = c36_cw(model,bank.S(linAll,:)); Sx = bank.S(linAll,:); jx = ev.Pj(linAll);
    pos = zeros(N*nt,1); pos(linAll) = 1:numel(linAll);
    v = sig(linAll).^2; Kps = cell(1,r); for l = 1:r, Kps{l} = zeros(numel(linAll),0); end
    selected = zeros(B,1);
    for b = 1:B
        s = sqrt(max(v,0)); rr = ev.G(linAll)./max(s,1e-300);
        Lcur = zeros(N,nt); Lcur(linAll) = s.^2.*0.5.*erfc(rr/sqrt(2)).*0.5.*erfc(-rr/sqrt(2));
        if b>1, Lcur(linAll(selected(1:b-1))) = 0; end
        R = inf(N,nt); R(linAll) = rr; rmin = min(R,[],2); rmin(~isfinite(rmin)) = 40;
        wTraj = (0.5*erfc(rmin/sqrt(2))).*(0.5*erfc(-rmin/sqrt(2)));      % Phi(-r^min)Phi(r^min)
        S2 = zeros(N,nt); S2(linAll) = s.^2;
        if b>1, S2(linAll(selected(1:b-1))) = 0; end
        [~,jt] = max(mean(S2.*wTraj,1)); [~,nb] = max(Lcur(:,jt));
        k = pos(nb + N*(jt-1));
        if k==0, [~,ia] = max(Lcur(linAll)); k = ia; end
        selected(b) = k;
        for l = 1:r
            cm = model.comp{l}; vstar = Vx(l,k);
            if vstar <= 1e-12*cm.sf2, continue; end
            kc = c36_ckern(cm,Sx,Sx(k,:)) - Wx{l}'*Wx{l}(:,k);
            if ~isempty(Kps{l}), kc = kc - Kps{l}*Kps{l}(k,:)'; end
            phi = kc/sqrt(vstar);
            Vx(l,:) = max(Vx(l,:) - (phi').^2,0);
            v = max(v - (model.Phi(jx,l).^2).*phi.^2,0);
            Kps{l} = [Kps{l}, phi];
        end
    end
    Snew = Sx(selected,:); info.selectedLin = linAll(selected); return;
end
if strcmp(method,'epak')
    lin = find(~isnan(ev.sig));
    rS = ev.G(lin)./max(ev.sig(lin),1e-12);         % signed ratio mu/sigma
    U = abs(rS);
    n = mod(lin-1,N)+1;
    rmin = accumarray(n,rS,[N 1],@min,Inf);          % trajectory-level minimum ratio
    % min-U point of every trajectory (snapshot location)
    best = Inf(N,1); bestLin = zeros(N,1);
    for q = 1:numel(lin)
        if U(q) < best(n(q)), best(n(q)) = U(q); bestLin(n(q)) = lin(q); end
    end
    tr = find(bestLin>0);
    pw = 0.5*erfc(abs(rmin(tr))/sqrt(2));            % p_w = Phi(-|r^min|) (EPAK)
    [~,o] = sort(pw,'descend'); tr = tr(o(1:min(Mtop,numel(tr)))); pw = pw(o(1:min(Mtop,numel(tr))));
    Sx = bank.S(bestLin(tr),:);
    % correlation of the dominant discrepancy component
    sfs = cellfun(@(x) x.sf2,model.comp); [~,ld] = max(sfs); cm = model.comp{ld};
    corrTo = @(A,Bm) c36_kernel(A,Bm,cm.ell,'m52');
    logIF = sum(log(max(1 - corrTo(Sx,model.S),1e-300)),2);
    selected = zeros(B,1);
    for b = 1:B
        sc = logIF + log(max(pw,1e-300));
        if b>1, sc(selected(1:b-1)) = -Inf; end
        [~,k] = max(sc); selected(b) = k;
        logIF = logIF + log(max(1 - corrTo(Sx,Sx(k,:)),1e-300));
    end
    Snew = Sx(selected,:); info.selectedLin = bestLin(tr(selected)); info.pw = pw(selected);
    return;
end
lin = find(~isnan(ev.sig));
mu = ev.G(lin); sg = max(ev.sig(lin),1e-12);
U = abs(mu)./sg;
switch method
    case {'ukb','pbalc','efflp'}
        cand = lin;
    case 'silk'
        failedRes = false(N,1);
        n = mod(lin-1,N)+1;
        hit = mu<0 & U>=2; failedRes(n(hit)) = true;
        keep = ~failedRes(n);
        cand = lin(keep); mu = mu(keep); sg = sg(keep); U = U(keep);
end
if isempty(cand)
    Snew = maximin(bank,model.S,B,rs); info.fallback = 'maximin'; return;
end
switch method
    case 'ukb',   score0 = -U;
    case 'silk',  score0 = -U;
    case 'pbalc', score0 = pbalcScore(mu,sg,cand,bank);
    case 'efflp', score0 = logEff(mu,sg);
end
[~,o] = sort(score0,'descend'); o = o(1:min(Mtop,numel(o)));
cand = cand(o); mu = mu(o); sg = sg(o);
Sx = bank.S(cand,:); jx = ev.Pj(cand);
selected = zeros(B,1); nsel = 0;
switch method
    case {'ukb','silk','pbalc'}
        % Kriging Believer: condition the variance on each selected block
        [Wx,zvx] = c36_cw(model,Sx); %#ok<ASGLU>
        v = sg.^2; Kps = cell(1,r); for l = 1:r, Kps{l} = zeros(numel(cand),0); end
        Vx = zvx;
        for b = 1:B
            s = sqrt(max(v,1e-300));
            switch method
                case {'ukb'}, sc = -abs(mu)./s;
                case 'silk'
                    sc = -abs(mu)./s;
                    % one candidate per trajectory: its min-U point
                    nn = mod(cand-1,N)+1; [~,~,g] = unique(nn);
                    best = accumarray(g,sc,[],@max); sc(sc < best(g)) = -Inf;
                case 'pbalc', sc = pbalcScore(mu,s,cand,bank);
            end
            if nsel>0, sc(ismember(Sx,Sx(selected(1:nsel),:),'rows')) = -Inf; end
            [~,k] = max(sc); nsel = nsel+1; selected(nsel) = k;
            for l = 1:r
                cm = model.comp{l}; vstar = Vx(l,k);
                if vstar <= 1e-12*cm.sf2, continue; end
                kc = c36_ckern(cm,Sx,Sx(k,:)) - Wx{l}'*Wx{l}(:,k);
                if ~isempty(Kps{l}), kc = kc - Kps{l}*Kps{l}(k,:)'; end
                phi = kc/sqrt(vstar);
                Vx(l,:) = max(Vx(l,:) - (phi').^2,0);
                v = max(v - (model.Phi(jx,l).^2).*phi.^2, 0);
                Kps{l} = [Kps{l}, phi];
            end
        end
    case 'efflp'
        % Lipschitz constant of the mean margin in normalised state space
        sc0 = [6 6 4];
        np = min(64,numel(cand)); ip = rs.randperm(numel(cand),np); probe = Sx(ip,:); jp = jx(ip);
        h = 1e-4; Lg = zeros(np,1); grad = zeros(np,3);
        for dd = 1:3
            e = zeros(1,3); e(dd) = h*sc0(dd);
            mp = c36_cgp_predict(model,probe+e); mm = c36_cgp_predict(model,probe-e);
            % gradient of the margin min_p g at the spatial minimiser (envelope theorem)
            idx = sub2ind(size(mp),jp(:)',1:np);
            grad(:,dd) = abs(mp(idx)-mm(idx))'/(2*h);
        end
        Lg(:) = sqrt(sum(grad.^2,2)); L = max(Lg);
        pen = zeros(numel(cand),1); Xn = Sx./sc0;
        for b = 1:B
            sc = logEff(mu,sg) + pen;
            if nsel>0, sc(ismember(Sx,Sx(selected(1:nsel),:),'rows')) = -Inf; end
            [~,k] = max(sc); nsel = nsel+1; selected(nsel) = k;
            rr = sqrt(sum((Xn - Xn(k,:)).^2,2));
            pen = pen + logPenalty(rr,mu(k),sg(k),L);
        end
        info.L = L;
end
Snew = Sx(selected(1:nsel),:);
info.selectedLin = cand(selected(1:nsel));
end

function s = pbalcScore(mu,sg,lin,bank)
% log{[Phi(b-U)-Phi(-b-U)] f(x)}, b = 1; f: density of the trajectory inputs
% (zeta, m) and of the instantaneous load state (Z1,Z2 | t).
b = 1; U = abs(mu)./sg;
pr = 0.5*erfc((U-b)/sqrt(2)) - 0.5*erfc((U+b)/sqrt(2));
N = bank.N; n = mod(lin-1,N)+1; S = bank.S(lin,:);
gam = 0.9 + 0.2*sin(pi*S(:,3)/4);
ld = -0.5*bank.zeta(n).^2 - 0.5*bank.um(n).^2 - 0.5*S(:,1).^2 - 0.5*(S(:,2)./gam).^2 - log(gam);
s = log(max(pr,1e-300)) + ld;
end

function v = logEff(mu,sg)
% Expected feasibility function (Bichon et al. 2008) with zbar = 0, eps = 2 sigma
e = 2*sg; t0 = -mu./sg; tm = (-e-mu)./sg; tp = (e-mu)./sg;
Ph = @(x) 0.5*erfc(-x/sqrt(2)); ph = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
eff = mu.*(2*Ph(t0) - Ph(tm) - Ph(tp)) - sg.*(2*ph(t0) - ph(tm) - ph(tp)) + e.*(Ph(tp) - Ph(tm));
v = log(max(eff,1e-300));
end

function v = logPenalty(r,muj,sgj,L)
% log P(-L r < G_j < L r), G_j ~ N(muj, sgj^2)
hi = (L*r - muj)/sgj; lo = (-L*r - muj)/sgj;
p = 0.5*erfc(-hi/sqrt(2)) - 0.5*erfc(-lo/sqrt(2));
v = log(max(p,1e-300));
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
