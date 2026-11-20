function [N,Nx,Nxx] = c2_hermite(r,le)
%C2_HERMITE Hermite-cubic shape functions at the element coordinate r in [0,1];
%   DOF = [wL thetaL wR thetaR]. The rotation is theta = dw/dx, so N2 and N4
%   contain the physical element length le.
r=r(:);
N=[1-3*r.^2+2*r.^3,le*(r-2*r.^2+r.^3), ...
    3*r.^2-2*r.^3,le*(-r.^2+r.^3)];
Nx=[(-6*r+6*r.^2)/le,1-4*r+3*r.^2, ...
    (6*r-6*r.^2)/le,-2*r+3*r.^2];
Nxx=[(-6+12*r)/le^2,(-4+6*r)/le, ...
    (6-12*r)/le^2,(-2+6*r)/le];
end
