function [Z1,Z2] = c36_process(xi,t)
%C36_PROCESS  Two Gaussian load processes from 8 standard-normal coefficients.
%   xi : N-by-8, one row per load history (shared by all times)
%   t  : 1-by-nt row of times
%   Z1, Z2 : N-by-nt
t = t(:)';
w = pi/2;
b = [cos(w*t); sin(w*t); 0.5*cos(2*w*t); 0.5*sin(2*w*t)]/sqrt(1.25);
Z1 = xi(:,1:4)*b;
Z2 = (xi(:,5:8)*b).*(0.9+0.2*sin(pi*t/4));
end
