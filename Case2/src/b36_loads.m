function ld = b36_loads(xi)
%B36_LOADS  Design-independent load quantities of histories (rows of xi).
%   OPERATOR LEVEL (exact): the beam is statically determinate, so the
%   quasi-static bending moment solves M'' = -q, M(0) = M(L) = 0, and does
%   not depend on the design, the section or the material:
%       Mqs(x,t) = int_0^L G(x,s) q(s,t) ds,  G = Green's function.
%   It is computed from the load history alone (no FE, no dynamics).
%   The modal forces of the assumed sine modes, f_j(t) = int q phi_j ds,
%   are design-independent as well.
%   ld.Mqs : nx-by-ntw-by-n (single, N m) on the retained window
%   ld.f   : n-by-nLF-by-J (single, N) on the LF integration grid c.tLF
persistent Gw Phi
c = b36_config(); cfg = c.hf; L = c.L;
if isempty(Gw)
    xs = c.xs; dx = diff(xs); ws = zeros(numel(xs),1);
    ws(1:end-1) = ws(1:end-1) + dx/2; ws(2:end) = ws(2:end) + dx/2;
    [X,S] = ndgrid(c.xw,xs');
    G = (S.*(L-X)/L).*(S<=X) + (X.*(L-S)/L).*(S>X);
    Gw = G.*ws';                                    % nx-by-ns
    Phi = sin(pi*xs*(1:c.J)/L).*ws;                  % ns-by-J
end
n = size(xi,1); nLF = numel(c.tLF);
Mqs = zeros(c.nx,c.ntw,n,'single'); f = zeros(n,nLF,c.J,'single');
for k = 1:n
    q = c2_load(c.xs,c.tLF,xi(k,:),cfg);            % ns-by-nLF
    Mqs(:,:,k) = single(Gw*q(:,c.twInLF));
    f(k,:,:) = reshape(single((Phi'*q)'),1,nLF,c.J);
end
ld = struct('Mqs',Mqs,'f',f,'n',n);
end
