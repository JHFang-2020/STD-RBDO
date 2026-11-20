function g = c36_phi2gain(a,rho)
%C36_PHI2GAIN  Expected reduction of the indicator variance pi(1-pi).
%   For a Gaussian margin M ~ N(mu,v), pi = Phi(a), a = -mu/sqrt(v), and a
%   prospective observation explaining a fraction rho = s^2/v of v:
%     E[pi(1-pi)] - E[pi+(1-pi+)] = Phi2(a,a;rho) - Phi(a)^2
%                                 = (1/2pi) int_0^{asin rho} exp(-a^2/(1+sin th)) dth
%   (Chevalier et al. 2014, Technometrics 56:455; exact for fixed GP params).
%   a : E-by-1 (or same size as rho); rho : E-by-n.  Gauss-Legendre (16 nodes).
persistent x w
if isempty(x)
    n = 16; j = (1:n-1)'; b = j./sqrt(4*j.^2-1);
    [V,D] = eig(diag(b,1)+diag(b,-1));
    [x,ix] = sort(diag(D)); w = 2*V(1,ix)'.^2;
    x = (x+1)/2; w = w/2;                 % nodes/weights on [0,1]
end
th = asin(min(max(rho,0),1-1e-12));
a2 = a.^2;
g = zeros(size(rho));
for q = 1:numel(x)
    g = g + w(q)*exp(-a2./(1+sin(th*x(q))));
end
g = g.*th/(2*pi);
end
