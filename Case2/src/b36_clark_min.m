function [m,v,w] = b36_clark_min(mu,C)
%b36_clark_min  Moment-matched Gaussian for min_k X_k, X ~ N(mu,C) (Clark 1961).
%   Returns mean m, variance v and tightness weights w (sum 1) such that
%   Cov(min X, Y) ~= sum_k w_k Cov(X_k, Y) for any Y jointly Gaussian with X
%   (exact for two variables).  Variables are combined in increasing mean.
mu = mu(:); K = numel(mu);
[~,ord] = sort(mu);
k1 = ord(1);
m = mu(k1); v = max(C(k1,k1),0);
w = zeros(K,1); w(k1) = 1;
cc = C(:,k1);                         % Cov(current min, X_k)
for ii = 2:K
    j = ord(ii); vj = max(C(j,j),0);
    th2 = v + vj - 2*cc(j);
    if th2 <= 1e-14*max(v+vj,realmin)
        if mu(j) < m, m = mu(j); v = vj; w(:) = 0; w(j) = 1; cc = C(:,j); end
        continue;
    end
    th = sqrt(th2);
    b = (mu(j) - m)/th;
    P = 0.5*erfc(-b/sqrt(2));          % P(current < X_j)
    ph = exp(-0.5*b^2)/sqrt(2*pi);
    m1 = m*P + mu(j)*(1-P) - th*ph;
    s2 = (m^2 + v)*P + (mu(j)^2 + vj)*(1-P) - (m + mu(j))*th*ph;
    v = max(s2 - m1^2,0); m = m1;
    cc = P*cc + (1-P)*C(:,j);
    w = P*w; w(j) = w(j) + (1-P);
end
end

