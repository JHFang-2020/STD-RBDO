function P = c36_pointfield(pop,bank,d)
%C36_POINTFIELD  Spatial minimum of the posterior-mean margin at every (n,t).
%   g(n,t,p) = R(d,zeta_n,p) - mu_D(s_nt,p)   (common mode NOT included)
%   P.g    : N-by-nt mean of min_p g (parabolic refinement in p)
%   P.j    : N-by-nt nearest retained node of the minimiser (for variances)
%   P.pst  : N-by-nt refined minimiser
c = c36_config(); a = c.a; p = c.pGrid'; np = numel(p);
N = bank.N; nt = bank.nt;
R = c.offset + (d(1) + a*bank.zeta).*(1-p) + (d(2) - a*bank.zeta).*p;   % N-by-np
F = reshape(R',np,N,1) - double(pop.MU);                                  % np-by-N-by-nt
[fmin,j] = min(F,[],1); fmin = reshape(fmin,N,nt); j = reshape(j,N,nt);
g = fmin; pst = p(j);
in = j>1 & j<np;
if any(in(:))
    [n1,t1] = find(in); jj = j(in);
    f0 = F(sub2ind([np N nt],jj-1,n1,t1)); f2 = F(sub2ind([np N nt],jj+1,n1,t1)); f1 = fmin(in);
    den = f0 - 2*f1 + f2; ok = den > 1e-14;
    delta = zeros(size(f1)); delta(ok) = 0.5*(f0(ok)-f2(ok))./den(ok);
    delta = max(min(delta,0.5),-0.5);
    g(in) = f1 - 0.25*(f0-f2).*delta;
    pst(in) = p(jj)' + delta*(p(2)-p(1));
end
P = struct('g',g,'j',j,'pst',pst,'d',d);
end
