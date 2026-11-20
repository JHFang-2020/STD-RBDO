function model = c36_cgp_fit(S,Y,spec,prev)
%C36_CGP_FIT  Component form of the spatiotemporal (physics-informed) Kriging.
%   Every snapshot returns a block Y(:,i) = D(s_i, pGrid).  The field is
%       D(s,p) = mu0(p) + sum_l Phi(p,l) z_l(s),   z_l ~ independent GPs
%   spec.spatial = 'op' : PHYSICS OPERATOR.  The quasi-static spatial
%                         operator L = d^4/dp^4 with L[D] = -96 (beam-type
%                         uniform-load equation) is imposed exactly: mu0 is
%                         the particular solution -4p^4 and Phi spans the
%                         null space {1,p,p^2,p^3} (orthonormalised on the
%                         grid).  Linearly-constrained GP (Jidling et al. 2017).
%                = 'se' : no operator: separable SE kernel in p (Kronecker),
%                         components = eigenvectors of K_p.
%   spec.mean    = 'lf' : calibrated low-fidelity physics trend + discrepancy
%                         z_l = b0_l + b1_l z_l^LF(s) + delta_l(s)  (KOH-type)
%                = 'const': constant trend (ordinary Kriging)
%                = 'lfonly': physics only, no learning.
%   S : Ns-by-3 states; Y : np-by-Ns blocks.  prev : previous model (warm start).
if nargin<4, prev = []; end
c = c36_config(); p = c.pGrid; np = numel(p);
model = struct('spec',spec,'S',S,'Y',Y,'pGrid',p);
switch spec.spatial
    case 'op'
        mu0 = -c.lambda*p.^4;
        [Phi,~] = qr([ones(np,1) p p.^2 p.^3],0);
    case 'se'
        mu0 = zeros(np,1); Phi = [];
    otherwise, error('c36_cgp_fit:spatial','unknown spatial model');
end
if strcmp(spec.mean,'lfonly')
    if isempty(Phi), Phi = eye(np); end
    model.mu0 = mu0; model.Phi = Phi; model.r = size(Phi,2); model.lfOnly = true; return;
end
model.lfOnly = false;
if strcmp(spec.spatial,'se')
    % shared Kronecker fit, then expressed as independent eigen-components
    if isempty(prev), pk = []; else, pk = prev.kron; end
    km = c36_gp_fit(S,Y,ternary(strcmp(spec.mean,'lf'),'pik','ok'),pk);
    keep = km.lp > 1e-10*km.lp(1);
    Phi = km.Up(:,keep); r = size(Phi,2);
    model.kron = km; model.mu0 = mu0; model.Phi = Phi; model.r = r;
    Zd = Phi'*Y;                                  % r-by-Ns data components
    if strcmp(spec.mean,'lf')
        T = km.beta(1)*Phi'*ones(np,1) + km.beta(2)*Phi'*c36_loadeffect(S,p,'lf');   % r-by-Ns
    else
        T = repmat(km.beta(1)*Phi'*ones(np,1),1,size(S,1));
    end
    comps = cell(1,r);
    for l = 1:r
        sf2 = km.sigma2*km.lp(l);
        comps{l} = gpcomp(S,Zd(l,:)',km.ellS,sf2,km.nu*km.sigma2,T(l,:)',[],[]);
        if strcmp(spec.mean,'lf')
            comps{l}.trendSE = [km.beta(1)*sum(Phi(:,l)), km.beta(2)];
        else
            comps{l}.trendSE = [km.beta(1)*sum(Phi(:,l)), 0];
        end
    end
    model.comp = comps;
    return;
end
% ---- operator model: independent component GPs on the null-space basis
r = size(Phi,2); model.mu0 = mu0; model.Phi = Phi; model.r = r;
Zd = Phi'*(Y - mu0);                              % r-by-Ns, exact projection
if any(strcmp(spec.mean,{'lf','lf1'}))
    Zlf = Phi'*(c36_loadeffect(S,p,'lf') - mu0);
end
comps = cell(1,r);
for l = 1:r
    off = zeros(size(S,1),1);
    switch spec.mean
        case 'lf',  H = [ones(size(S,1),1), Zlf(l,:)'];   % calibrated scale + offset
        case 'lf1', H = ones(size(S,1),1); off = Zlf(l,:)'; % residual on D - D_LF (manuscript Eq. 4)
        otherwise,  H = ones(size(S,1),1);
    end
    th0 = []; if ~isempty(prev), th0 = prev.comp{l}.theta; end
    ns = isfield(spec,'nonstat') && spec.nonstat;
    comps{l} = fitcomp(S,Zd(l,:)' - off,H,th0,ns);
    comps{l}.offsetLF = any(off~=0) || strcmp(spec.mean,'lf1');
end
model.comp = comps;
end

function cmp = fitcomp(S,z,H,th0,ns)
% theta = [log ell (3)] or, for the non-stationary amplitude, [log ell, log gamma]
lb = log([0.25 0.25 0.15]); ub = log([40 40 30]);
starts = [log([2 2 1.5]); log([1 1 0.8]); log([6 6 3])];
if ns
    lb = [lb log(1e-3)]; ub = [ub log(10)];
    starts = [starts, log([0.1; 0.5; 0.02])];
end
if ~isempty(th0)
    if numel(th0) < numel(lb), th0 = [th0 log(0.1)]; end
    if numel(th0) > numel(lb), th0 = th0(1:numel(lb)); end
    starts = [th0; starts(1,:)];
end
fo = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',300,'OptimalityTolerance',1e-6,'StepTolerance',1e-7);
best = Inf; th = starts(1,:); nuRel = 1e-8;
for k = 1:size(starts,1)
    x0 = min(max(starts(k,:),lb+1e-6),ub-1e-6);
    try
        [x,f] = fmincon(@(t) cnll(t,S,z,H,nuRel),x0,[],[],[],[],lb,ub,[],fo);
    catch
        x = x0; f = cnll(x0,S,z,H,nuRel);
    end
    if f < best, best = f; th = x; end
end
[~,sf2,beta] = cnll(th,S,z,H,nuRel);
gam = 0; if numel(th) == 4, gam = exp(th(4)); end
cmp = gpcomp(S,z,exp(th(1:3)),sf2,nuRel*sf2,H*beta,H,beta,gam);
cmp.theta = th; cmp.nll = best;
end

function [f,sf2,beta] = cnll(th,S,z,H,nuRel)
n = numel(z);
K = c36_kernel(S,S,exp(th(1:3)),'m52');
if numel(th) == 4
    w = c36_wamp(S,exp(th(4))); K = K.*(w*w');
end
K = K + nuRel*eye(n);
[L,pd] = chol(K,'lower');
if pd, f = 1e20; sf2 = NaN; beta = zeros(size(H,2),1); return; end
Hi = L\H; zi = L\z;
beta = (Hi'*Hi + 1e-12*eye(size(H,2)))\(Hi'*zi);
r = zi - Hi*beta;
sf2 = max(r'*r/n,1e-300);
f = 0.5*n*log(sf2) + sum(log(diag(L)));
end

function cmp = gpcomp(S,z,ell,sf2,nugget,trendTrain,H,beta,gam)
if nargin<9, gam = 0; end
n = size(S,1);
cmp = struct('S',S,'ell',ell,'sf2',sf2,'gamma',gam);
K = c36_ckern(cmp,S,S) + nugget*eye(n);
L = chol((K+K')/2,'lower');
cmp.nugget = nugget; cmp.L = L; cmp.alpha = L'\(L\(z - trendTrain));
cmp.H = H; cmp.beta = beta; cmp.theta = log(ell);
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
