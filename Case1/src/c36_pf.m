function [pf,grad,se,q] = c36_pf(E,pGrid,d,zeta)
%C36_PF  Trajectory failure probability with the common mode integrated exactly.
%   Pf(d) = mean_n Phi(-h_n(d)/a)
%   grad  : 1-by-2 analytic gradient;  se : sample s.d./sqrt(N) of q_n
%   q     : conditional failure probability of every trajectory.
c = c36_config(); a = c.a;
[h,~,dh] = c36_margin(E,pGrid,d,zeta);
q = 0.5*erfc(h/(sqrt(2)*a));
pf = mean(q);
if nargout>1
    dens = exp(-0.5*(h/a).^2)/(sqrt(2*pi)*a);
    grad = -mean(dens.*dh,1);
end
if nargout>2
    se = std(q)/sqrt(numel(q));
end
end
