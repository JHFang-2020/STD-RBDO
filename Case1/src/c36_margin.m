function [h,pstar,dh] = c36_margin(E,pGrid,d,zeta)
%C36_MARGIN  Trajectory margin without the common mode, from a time envelope.
%   h_n(d) = min_p [0.35 + (d1+a*zeta_n)(1-p) + (d2-a*zeta_n)p - e_n(p)]
%   so that  g_min,n = m_n + h_n(d),  m_n ~ N(0,a^2)  (exact).
%   The discrete minimum is refined by a three-point parabola.
%   dh : N-by-2 gradient of h with respect to d (envelope theorem).
c = c36_config(); a = c.a;
p = pGrid(:)'; np = numel(p); N = size(E,1);
R = c.offset + (d(1) + a*zeta(:)).*(1-p) + (d(2) - a*zeta(:)).*p;   % N-by-np
F = R - E;
[fmin,j] = min(F,[],2);
h = fmin; pstar = p(j)';
in = j>1 & j<np;
if any(in)
    jj = j(in); rows = find(in);
    f0 = F(sub2ind([N np],rows,jj-1));
    f1 = fmin(in);
    f2 = F(sub2ind([N np],rows,jj+1));
    den = f0 - 2*f1 + f2;
    ok = den > 1e-14;
    delta = zeros(size(f1));
    delta(ok) = 0.5*(f0(ok)-f2(ok))./den(ok);
    delta = max(min(delta,0.5),-0.5);
    hp = p(2)-p(1);
    h(in)     = f1 - 0.25*(f0-f2).*delta;
    pstar(in) = p(jj)' + delta*hp;
end
if nargout>2
    dh = [1-pstar, pstar];
end
end
