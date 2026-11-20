function [E,tStar] = c36_true_envelope(xi,pGrid,tGrid,chunk)
%C36_TRUE_ENVELOPE  Exact time envelope e_n(p) = max_t D_hf(s_n(t),p).
%   Because the capacity R(X,p) does not depend on t,
%       min_{p,t} G = min_p [ R(X,p) - e_n(p) ],
%   so the envelope on a fine (p,t) grid is a design-free sufficient statistic
%   of one load history.  Used ONLY for benchmark/validation, never for learning.
%   xi : N-by-8;  E : N-by-numel(pGrid);  tStar : time index of the max (N-by-np)
if nargin<4, chunk = 256; end
N = size(xi,1); np = numel(pGrid); nt = numel(tGrid);
E = zeros(N,np); tStar = zeros(N,np,'uint16');
for s = 1:chunk:N
    ii = s:min(N,s+chunk-1); n = numel(ii);
    [Z1,Z2] = c36_process(xi(ii,:),tGrid);
    S = [Z1(:), Z2(:), repelem(tGrid(:),n)];
    D = c36_loadeffect(S,pGrid,'hf');            % np-by-(n*nt), time-major
    D = reshape(D,np,n,nt);
    [m,it] = max(D,[],3);
    E(ii,:) = m'; tStar(ii,:) = uint16(it');
end
end
