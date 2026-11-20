function [mu,v,zm,zv] = c36_cgp_predict(model,Ss)
%C36_CGP_PREDICT  Posterior of the load-effect blocks at states Ss (M-by-3).
%   mu, v : np-by-M mean and variance on the retained spatial grid
%   zm, zv: r-by-M component means/variances (independent components)
p = model.pGrid; M = size(Ss,1); r = model.r;
if model.lfOnly
    mu = c36_loadeffect(Ss,p,'lf'); v = zeros(size(mu));
    zm = model.Phi'*(mu - model.mu0); zv = zeros(size(zm)); return;
end
needLF = any(strcmp(model.spec.mean,{'lf','lf1'}));
if needLF
    Dlf = c36_loadeffect(Ss,p,'lf');
    if strcmp(model.spec.spatial,'op'), Zlf = model.Phi'*(Dlf - model.mu0);
    else, Zlf = model.Phi'*Dlf; end
end
zm = zeros(r,M);
wantVar = nargout>1;
if wantVar, zv = zeros(r,M); end
for l = 1:r
    cm = model.comp{l};
    if strcmp(model.spec.spatial,'op')
        switch model.spec.mean
            case 'lf',  tr = cm.beta(1) + cm.beta(2)*Zlf(l,:);
            case 'lf1', tr = cm.beta(1) + Zlf(l,:);
            otherwise,  tr = cm.beta(1)*ones(1,M);
        end
    else
        tr = cm.trendSE(1)*ones(1,M);
        if needLF, tr = tr + cm.trendSE(2)*Zlf(l,:); end
    end
    Kxs = c36_ckern(cm,cm.S,Ss);                        % Ns-by-M
    zm(l,:) = tr + cm.alpha'*Kxs;
    if wantVar
        W = cm.L\Kxs;
        zv(l,:) = max(priorvar(cm,Ss) - sum(W.^2,1),0);
    end
end
mu = model.mu0 + model.Phi*zm;
if wantVar, v = (model.Phi.^2)*zv; end
end

function pv = priorvar(cm,Ss)
pv = cm.sf2*ones(1,size(Ss,1));
if isfield(cm,'gamma') && cm.gamma > 0, pv = pv.*(c36_wamp(Ss,cm.gamma).^2)'; end
end
