function [Snew,info] = c36_acq_rv(model,bank,funs,dw,B,opts)
%C36_ACQ_RV  ST-MSUR-RV: batch SUR on the posterior variance of the working risk.
%   funs : cell of linearised risk functionals (c36_rvfun) at the designs of
%          the neighbourhood, weights dw.  A snapshot at state c observes the
%          block z_l(c), l=1..r; the variance of L_j it explains is
%          sum_l Cov(L_j,z_l(c))^2 / Var(z_l(c)) (independent components),
%          and later batch members are scored conditionally on the earlier
%          ones through exact Schur-complement updates (variances only).
%   Criterion:  max_c  sum_j dw_j * explained_j(c) / Var(L_j)
if nargin<6, opts = struct(); end
nExplore = getf(opts,'nExplore',200);
rs = RandStream('mt19937ar','Seed',getf(opts,'seed',1));
r = model.r; J = numel(funs);
candLin = [];
for j = 1:J, candLin = [candLin; funs{j}.lin]; end %#ok<AGROW>
expl = rs.randperm(bank.N*bank.nt,min(nExplore,bank.N*bank.nt))';
candLin = unique([candLin; expl]);
Sc = bank.S(candLin,:); nC = numel(candLin);
[Wc,Vc] = c36_cw(model,Sc);
XL = cell(1,r); for l = 1:r, XL{l} = zeros(J,nC); end
varL = zeros(J,1);
for j = 1:J
    f = funs{j}; varL(j) = max(f.varL,1e-300);
    [Wf,~] = c36_cw(model,f.S);
    for l = 1:r
        cm = model.comp{l}; al = f.coef.*model.Phi(f.j,l);
        kv = zeros(1,nC);
        for s = 1:2000:numel(al)
            ii = s:min(numel(al),s+1999);
            kv = kv + al(ii)'*c36_ckern(cm,f.S(ii,:),Sc);
        end
        XL{l}(j,:) = kv - (Wf{l}*al)'*Wc{l};
    end
end
wj = dw(:)./varL;
selected = zeros(B,1); scores = zeros(B,1); Kcs = cell(1,r);
for l = 1:r, Kcs{l} = zeros(nC,0); end
explained = zeros(J,1);
for b = 1:B
    sc = zeros(1,nC);
    for l = 1:r
        ok = Vc(l,:) > 1e-12*model.comp{l}.sf2;
        sc(ok) = sc(ok) + (wj'*(XL{l}(:,ok).^2))./Vc(l,ok);
    end
    if b>1, sc(selected(1:b-1)) = -Inf; end
    [scores(b),cs] = max(sc); selected(b) = cs;
    for l = 1:r
        cm = model.comp{l}; vstar = Vc(l,cs);
        if vstar <= 1e-12*cm.sf2, continue; end
        kc = c36_ckern(cm,Sc,Sc(cs,:)) - Wc{l}'*Wc{l}(:,cs);
        if ~isempty(Kcs{l}), kc = kc - Kcs{l}*Kcs{l}(cs,:)'; end
        phiC = kc/sqrt(vstar); phiL = XL{l}(:,cs)/sqrt(vstar);
        explained = explained + phiL.^2;
        Vc(l,:) = max(Vc(l,:) - (phiC').^2,0);
        XL{l} = XL{l} - phiL*phiC';
        Kcs{l} = [Kcs{l}, phiC];
    end
end
Snew = Sc(selected,:);
info = struct('scores',scores,'selectedLin',candLin(selected),'nCandidates',nC, ...
    'varL',varL,'explainedFraction',explained./varL);
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
