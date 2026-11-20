function gp = b36_gp_fit(Z,y,spec,prev)
%B36_GP_FIT  Physics-informed discrepancy GP of the beam (Case 2, v36).
%   y = log W_HF - log W_mean at critical windows (one row per window),
%   W_mean = simplified-physics window maximum (spec.mean = 'lf') or the
%   exact quasi-static one ('qs').  Prior:
%       delta(z) ~ GP( beta, sf2 * w(z) w(z') * k_M52(z,z') ) + nugget,
%       w(z)^2 = 1 + gamma * a(z)^2   (non-stationary amplitude, spec.nonstat)
%   with grouped ARD length scales [d1 d2 d3 | xi_manuf | xi_mat | xi_load | t | a],
%   beta by GLS, sf2 profiled, (ell, gamma, nugget) by maximum likelihood.
%   Z : n-by-22 (b36_feat).  prev : previous gp (warm start).
if nargin<4, prev = []; end
grp = [1 2 3, 4 4 4, 5 5, 6*ones(1,12), 7, 8];
ng = 8; n = size(Z,1);
ns = isfield(spec,'nonstat') && spec.nonstat;
D2 = zeros(n,n,ng);
for g = 1:ng
    A = Z(:,grp==g); D2(:,:,g) = max(sum(A.^2,2) + sum(A.^2,2)' - 2*(A*A'),0);
end
lb = log([0.02 0.02 0.02 0.3 0.3 0.3 0.01 0.01]); ub = log([20 20 20 50 50 50 10 10]);
th0 = log([0.3 0.3 0.3 3 3 3 0.05 0.3]);
lb = [lb log(1e-6)]; ub = [ub log(0.5)]; th0 = [th0 log(1e-3)];            % nugget (relative)
if ns, lb = [lb log(1e-3)]; ub = [ub log(1e3)]; th0 = [th0 log(1)]; end    % gamma
starts = th0;
if ~isempty(prev) && numel(prev.theta)==numel(th0), starts = [prev.theta; th0]; end
fo = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',400,'OptimalityTolerance',1e-5,'StepTolerance',1e-6);
best = Inf; th = starts(1,:);
for k = 1:size(starts,1)
    x0 = min(max(starts(k,:),lb+1e-6),ub-1e-6);
    try
        [x,f] = fmincon(@(t) nll(t,D2,Z(:,end),y,ns),x0,[],[],[],[],lb,ub,[],fo);
    catch
        x = x0; f = nll(x0,D2,Z(:,end),y,ns);
    end
    if f < best, best = f; th = x; end
end
[~,sf2,beta,L,alpha] = nll(th,D2,Z(:,end),y,ns);
gam = 0; if ns, gam = exp(th(10)); end
gp = struct('Z',Z,'y',y,'grp',grp,'ell',exp(th(1:8)),'nu',exp(th(9)),'gamma',gam, ...
    'sf2',sf2,'beta',beta,'L',L,'alpha',alpha,'theta',th,'nll',best,'spec',spec,'n',n);
gp.tau2 = gp.nu*sf2;                                  % observation (nugget) variance
end

function [f,sf2,beta,L,alpha] = nll(th,D2,a,y,ns)
n = numel(y);
R2 = zeros(n);
for g = 1:size(D2,3), R2 = R2 + D2(:,:,g)/exp(2*th(g)); end
r = sqrt(5*R2); K = (1 + r + r.^2/3).*exp(-r);
if ns
    w = sqrt(1 + exp(th(10))*a.^2); K = K.*(w*w');
end
K = K + exp(th(9))*eye(n);
[L,pd] = chol(K,'lower');
if pd, f = 1e20; sf2 = NaN; beta = 0; alpha = []; return; end
o = ones(n,1); Lo = L\o; Ly = L\y;
beta = (Lo'*Ly)/(Lo'*Lo);
res = Ly - Lo*beta;
sf2 = max(res'*res/n,1e-300);
f = 0.5*n*log(sf2) + sum(log(diag(L)));
if nargout>3
    L = sqrt(sf2)*L;                                  % chol of sf2*(K+nu I)
    alpha = L'\(L\(y - beta));
end
end
