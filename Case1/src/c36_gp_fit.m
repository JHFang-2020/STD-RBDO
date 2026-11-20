function model = c36_gp_fit(S,Y,meanType,prev,opts)
%C36_GP_FIT  Physics-informed spatiotemporal Kriging for snapshot blocks.
%   D(s,p) = beta' h(s,p) + delta(s,p),   delta ~ GP(0, sigma^2 k_s(s,s') k_p(p,p'))
%   meanType 'pik' : h = [1, D_LF(s,p)]   (physics mean with calibrated
%                    scale/offset + discrepancy; KOH-type linear calibration)
%            'ok'  : h = 1                (ordinary Kriging, no physics)
%            'lf'  : D = D_LF             (physics only, no learning)
%   Every snapshot returns the full retained block on the common grid pGrid,
%   so the data covariance is sigma^2 (K_s kron K_p + nu I).  Exact inference
%   uses the eigen-decompositions of the two factors (Saatci 2012; Stegle et
%   al. 2011): cost O(Ns^3 + np^3) instead of O((Ns*np)^3).
%   S : Ns-by-3 states, Y : np-by-Ns observed blocks.
if nargin<4, prev = []; end
if nargin<5, opts = struct(); end
c = c36_config(); p = c.pGrid;
model = struct('meanType',meanType,'S',S,'Y',Y,'pGrid',p,'nu',1e-8);
if strcmp(meanType,'lf')
    model.lfOnly = true; return;
end
model.lfOnly = false;
switch meanType
    case 'pik', H = {ones(size(Y)), c36_loadeffect(S,p,'lf')};
    case 'ok',  H = {ones(size(Y))};
    otherwise, error('c36_gp_fit:mean','unknown mean type');
end
lb = log([0.25 0.25 0.15 0.04]); ub = log([30 30 20 3]);
if isempty(prev)
    starts = [log([2 2 1.5 0.3]); log([1 1 0.8 0.2]); log([5 5 3 0.5]); log([0.7 1.5 0.5 0.25])];
else
    starts = [prev.theta; log([2 2 1.5 0.3])];
end
if isfield(opts,'nStarts'), starts = starts(1:min(end,opts.nStarts),:); end
fo = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',400,'OptimalityTolerance',1e-6,'StepTolerance',1e-7);
best = Inf; theta = starts(1,:);
for k = 1:size(starts,1)
    x0 = min(max(starts(k,:),lb+1e-6),ub-1e-6);
    try
        [x,f] = fmincon(@(th) nll(th,S,Y,H,p,model.nu),x0,[],[],[],[],lb,ub,[],fo);
    catch
        x = x0; f = nll(x0,S,Y,H,p,model.nu);
    end
    if f < best, best = f; theta = x; end
end
[~,parts] = nll(theta,S,Y,H,p,model.nu);
model.theta = theta;
model.ellS = exp(theta(1:3)); model.ellP = exp(theta(4));
model.beta = parts.beta; model.sigma2 = parts.sigma2;
model.Up = parts.Up; model.lp = parts.lp; model.Us = parts.Us; model.ls = parts.ls;
model.Lam = parts.Lam; model.Kp = parts.Kp;
model.alpha = parts.Up*(parts.Rt./parts.Lam)*parts.Us';      % np-by-Ns
model.Malpha = parts.Kp*model.alpha;                          % mean weights
model.B = parts.lp.*parts.Up';                                % B(l,j) = lp_l Up(j,l)
model.invLam = 1./parts.Lam;
model.nll = best;
end

function [f,parts] = nll(th,S,Y,H,p,nu)
ellS = exp(th(1:3)); ellP = exp(th(4));
Kp = c36_kernel(p,p,ellP,'se'); Kp = (Kp+Kp')/2;
Ks = c36_kernel(S,S,ellS,'m52'); Ks = (Ks+Ks')/2;
[Up,lp] = eig(Kp,'vector'); [lp,ix] = sort(max(lp,0),'descend'); Up = Up(:,ix);
[Us,ls] = eig(Ks,'vector'); [ls,ix] = sort(max(ls,0),'descend'); Us = Us(:,ix);
Lam = lp*ls' + nu;
Yt = Up'*Y*Us;
nb = numel(H); Ht = cell(1,nb);
for k = 1:nb, Ht{k} = Up'*H{k}*Us; end
G = zeros(nb); b = zeros(nb,1);
for k = 1:nb
    b(k) = sum(Ht{k}.*Yt./Lam,'all');
    for l = k:nb
        G(k,l) = sum(Ht{k}.*Ht{l}./Lam,'all'); G(l,k) = G(k,l);
    end
end
beta = (G + 1e-12*trace(G)*eye(nb))\b;
Rt = Yt;
for k = 1:nb, Rt = Rt - beta(k)*Ht{k}; end
N = numel(Y);
sigma2 = max(sum(Rt.^2./Lam,'all')/N, 1e-300);
f = 0.5*N*log(sigma2) + 0.5*sum(log(Lam),'all');
if ~isfinite(f), f = 1e20; end
if nargout>1
    parts = struct('beta',beta,'sigma2',sigma2,'Up',Up,'lp',lp,'Us',Us,'ls',ls,'Lam',Lam,'Rt',Rt,'Kp',Kp);
end
end
