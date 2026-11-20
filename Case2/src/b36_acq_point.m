function [sel,info] = b36_acq_point(gp,E,P,B,method,seed,excl0)
%B36_ACQ_POINT  Literature batch learning functions on the shared surrogate.
%   Candidates: the active histories at the incumbent design (one HF call =
%   one history, observing its block of Kobs windows).  Learning functions
%   act on the tracked critical windows (n,k) with window margin
%   g = log(160/W_LF) - delta, posterior mean m and latent s.d. s.
%   'ukb'    : U = |m|/s (Echard et al. 2011) + Kriging Believer (Ginsbourger
%              et al. 2010: variance conditioned on the selected blocks, mean frozen)
%   'silk'   : SILK (Hu & Mahadevan 2016): trajectories with a window m<0, U>=2
%              are resolved as failed; value min_k U; Kriging Believer
%   'slabalq': SL-ABALQ (Dang et al. 2025, Eq. 32): time step (window) j* =
%              argmax_j mean_n s_nj^2 Phi(-r_n^min)Phi(r_n^min); sample n* =
%              argmax s^2 Phi(-r)Phi(r) at j*; Kriging Believer
%   'epak'   : EPAK (Hu et al. 2025): p_w = Phi(-|r^min|), batch by influence
%              function prod(1-R) over training and in-batch points times p_w
%              (posterior frozen within the batch: no Kriging Believer)
%   'pbalc'  : PBALC learning function 3 (Dang et al. 2024): [Phi(1-U)-Phi(-1-U)] f(xi);
%              Kriging Believer
%   'efflp'  : EFF (Bichon et al. 2008, eps = 2s) with local penalisation
%              P(|G_j| < L r) (Wang et al. 2025), L = largest difference quotient of the
%              window-margin mean m between nearest tracked windows; posterior frozen
%   excl0    : logical over E.idx, histories that may not be selected (exact repeats)
%   'random' : uniformly random active history (control)
%   E : b36_popeval at the incumbent (keepAll); returns sel = rows of E.idx
rs = RandStream('mt19937ar','Seed',seed);
n = numel(E.idx); K = E.K; Ko = b36_config().Kobs;
info = struct('method',method);
if n == 0, sel = []; return; end
if strcmp(method,'random')
    ok = (1:n)'; if nargin>=7 && ~isempty(excl0), ok = find(~excl0(:)); end
    sel = ok(rs.randperm(numel(ok),min(B,numel(ok)))); return;
end
m = E.m;
vt = zeros(n,K); for k = 1:K, vt(:,k) = squeeze(E.C(k,k,:)); end
C = b36_cands(E,P,(1:n)');
Zall = [E.Z; C.Z]; nt = n*K;                       % tracked windows, then candidate blocks
[~,~,Wall] = b36_gp_pred(gp,Zall);
Uall = zeros(size(Zall,1),0);
sel = zeros(B,1); excl = false(n,1); if nargin>=7 && ~isempty(excl0), excl = logical(excl0(:)); end
if strcmp(method,'efflp')
    Lg = lipschitz(gp,E.Z,m,rs);
    pen = zeros(n,K); info.L = Lg;
end
if strcmp(method,'epak')
    % snapshot location of each trajectory: its min-U window (as in EPAK)
    [~,kU] = min(abs(m)./sqrt(max(vt,1e-300)),[],2);
    Zp = E.Z((0:n-1)'*K + kU,:);
    logIF = logInfluence(gp,Zp,gp.Z);
end
for b = 1:B
    s = sqrt(max(vt,1e-300)); U = abs(m)./s; r = m./s;
    switch method
        case 'ukb'
            sc = -min(U,[],2);
        case 'silk'
            failed = any(m<0 & U>=2,2);
            sc = -min(U,[],2); sc(failed) = -Inf;
        case 'slabalq'
            rmin = min(r,[],2);
            wTraj = Phi(-rmin).*Phi(rmin);
            nW = b36_config().nW; S2 = zeros(n,nW);
            for k = 1:K
                S2(sub2ind([n nW],(1:n)',E.wk(:,k))) = s(:,k).^2;
            end
            S2(excl,:) = 0;
            [~,jt] = max(mean(S2.*wTraj,1));
            Lc = -Inf(n,1);
            for k = 1:K
                on = E.wk(:,k)==jt;
                Lc(on) = s(on,k).^2.*Phi(-r(on,k)).*Phi(r(on,k));
            end
            if all(~isfinite(Lc(~excl)))
                Lw = s.^2.*Phi(-r).*Phi(r); Lc = max(Lw,[],2);
            end
            sc = Lc;
        case 'epak'
            pw = Phi(-abs(min(r,[],2)));
            sc = logIF + log(max(pw,1e-300));
        case 'pbalc'
            lf = -0.5*sum(P.xi(E.idx,:).^2,2);
            pr = Phi(1-U) - Phi(-1-U);
            sc = max(log(max(pr,1e-300)),[],2) + lf;
        case 'efflp'
            sc = max(logEff(m,s) + pen,[],2);
        otherwise
            error('b36_acq_point:method','unknown method %s',method);
    end
    sc(excl) = -Inf;
    [best,i] = max(sc);
    if ~isfinite(best), sel = sel(1:b-1); break; end
    sel(b) = i; excl(i) = true;
    if ~any(strcmp(method,{'epak','efflp'}))      % EPAK / EFF+LP: batch diversity by IF / penalty only
        % ---- Kriging Believer: condition the variances on the selected block
        rowsY = nt + (i-1)*Ko + (1:Ko);
        Kc = b36_gp_k(gp,Zall,Zall(rowsY,:)) - Wall'*Wall(:,rowsY);
        if ~isempty(Uall), Kc = Kc - Uall*Uall(rowsY,:)'; end
        S = Kc(rowsY,:) + gp.tau2*eye(Ko); S = (S+S')/2;
        Ls = chol(S + 1e-12*mean(diag(S))*eye(Ko),'lower');
        Un = Kc/Ls'; Uall = [Uall, Un]; %#ok<AGROW>
        vt = max(vt - reshape(sum(Un(1:nt,:).^2,2),K,n)',0);
    end
    if strcmp(method,'epak')
        logIF = logIF + logInfluence(gp,Zp,C.Z((i-1)*Ko+(1:Ko),:));
    end
    if strcmp(method,'efflp')
        [~,kk] = max(logEff(m(i,:),s(i,:)));
        zj = E.Z((i-1)*K+kk,:); rr = normdist(gp,E.Z,zj);
        pen = pen + reshape(logPenalty(rr,m(i,kk),s(i,kk),Lg),K,n)';
    end
end
info.vtFinal = vt;                  % KB-conditioned window variances (unit test)
end

function p = Phi(x)
p = 0.5*erfc(-x/sqrt(2));
end

function li = logInfluence(gp,Zp,Zr)
% EPAK influence function: sum log(1-R) over reference points
% (R = prior correlation of the discrepancy GP)
k0 = b36_gp_kpair(gp,Zp,Zp);
kr = b36_gp_kpair(gp,Zr,Zr);
R = b36_gp_k(gp,Zp,Zr)./sqrt(k0*kr');
li = sum(log(max(1-R,1e-300)),2);
end

function r = normdist(gp,A,z)
r2 = zeros(size(A,1),1);
for g = 1:numel(gp.ell), r2 = r2 + sum(((A(:,gp.grp==g) - z(gp.grp==g))/gp.ell(g)).^2,2); end
r = sqrt(r2);
end

function L = lipschitz(gp,Z,m,rs)
% Lipschitz constant of the window-margin mean m = log(160/W_LF) - mu_delta in
% length-scale units: largest difference quotient between each sampled
% tracked window and its nearest other tracked window (the margin is not an
% explicit function of the features because W_LF comes from the physics)
mv = reshape(m',[],1);                        % row order of Z (history-major)
np = min(400,size(Z,1)); ip = rs.randperm(size(Z,1),np);
q = zeros(np,1);
for k = 1:np
    r = normdist(gp,Z,Z(ip(k),:)); r(ip(k)) = Inf;
    [rmin,j] = min(r);
    q(k) = abs(mv(ip(k)) - mv(j))/max(rmin,1e-12);
end
L = max(q);
end

function v = logEff(mu,sg)
% Expected feasibility function (Bichon et al. 2008), zbar = 0, eps = 2 sigma
e = 2*sg; t0 = -mu./sg; tm = (-e-mu)./sg; tp = (e-mu)./sg;
Ph = @(x) 0.5*erfc(-x/sqrt(2)); ph = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
eff = mu.*(2*Ph(t0) - Ph(tm) - Ph(tp)) - sg.*(2*ph(t0) - ph(tm) - ph(tp)) + e.*(Ph(tp) - Ph(tm));
v = log(max(eff,1e-300));
end

function v = logPenalty(r,muj,sgj,L)
hi = (L*r - muj)/sgj; lo = (-L*r - muj)/sgj;
p = 0.5*erfc(-hi/sqrt(2)) - 0.5*erfc(-lo/sqrt(2));
v = log(max(p,1e-300));
end
