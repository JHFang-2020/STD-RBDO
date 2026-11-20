function [W,zv] = c36_cw(model,Sx)
%C36_CW  Whitened cross-kernels of every component GP at states Sx.
%   W{l} = L_l \ k_l(S,Sx)   (Ns-by-nx), zv(l,:) posterior variances.
%   Posterior covariance of component l:  c_l(a,b) = k_l(a,b) - W{l}(:,a)'W{l}(:,b)
%   (k_l includes sf2 and the optional non-stationary amplitude, c36_ckern)
r = model.r; nx = size(Sx,1);
W = cell(1,r); zv = zeros(r,nx);
if model.lfOnly
    for l = 1:r, W{l} = zeros(0,nx); end
    return;
end
for l = 1:r
    cm = model.comp{l};
    W{l} = cm.L\c36_ckern(cm,cm.S,Sx);
    pv = cm.sf2*ones(1,nx);
    if isfield(cm,'gamma') && cm.gamma > 0, pv = pv.*(c36_wamp(Sx,cm.gamma).^2)'; end
    zv(l,:) = max(pv - sum(W{l}.^2,1),0);
end
end
