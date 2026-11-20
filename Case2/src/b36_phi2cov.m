function g = b36_phi2cov(a,b,rho)
%B36_PHI2COV  Phi2(a,b;rho) - Phi(a)Phi(b) (covariance of two Gaussian
%   threshold indicators), Gauss-Legendre (16 nodes) on
%   (1/2pi) int_0^{asin rho} exp(-(a^2+b^2-2ab sin t)/(2 cos^2 t)) dt.
%   a, b, rho : arrays of equal size (or scalar expansion).
persistent x w
if isempty(x)
    n = 16; j = (1:n-1)'; bb = j./sqrt(4*j.^2-1);
    [V,D] = eig(diag(bb,1)+diag(bb,-1));
    [x,ix] = sort(diag(D)); w = 2*V(1,ix)'.^2;
    x = (x+1)/2; w = w/2;
end
th = asin(min(max(rho,-1+1e-12),1-1e-12));
g = zeros(size(th));
for q = 1:numel(x)
    s = sin(th*x(q));
    g = g + w(q)*exp(-(a.^2 + b.^2 - 2*a.*b.*s)./(2*(1 - s.^2)));
end
g = g.*th/(2*pi);
end
