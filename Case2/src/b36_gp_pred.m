function [mu,v,W] = b36_gp_pred(gp,A)
%B36_GP_PRED  Posterior mean and latent variance of the discrepancy at rows A.
%   W = L \ k(Z,A) (whitened cross-covariance), so that the posterior
%   covariance of any two sets is  k(A,B) - W_A' W_B  (b36_gp_k).
Ks = b36_gp_k(gp,gp.Z,A);
mu = gp.beta + Ks'*gp.alpha;
if nargout>1
    W = gp.L\Ks;
    pv = gp.sf2*ones(size(A,1),1);
    if gp.gamma > 0, pv = pv.*(1 + gp.gamma*A(:,end).^2); end
    v = max(pv - sum(W.^2,1)',0);
end
end
