function Z = b36_feat(d,xi,tw,a)
%B36_FEAT  GP input of a critical window (one row per window).
%   [ (d-0.16)/0.14 (3) | xi (17) | (t_window-1)/3 | a ]
%   a = log(W_LF/W_QS): dynamic amplification of the window predicted by the
%   simplified physics (the LF is exact in the quasi-static limit a -> 0).
%   d : 1-by-3 or n-by-3; xi : n-by-17; tw, a : n-by-1
n = size(xi,1);
if size(d,1)==1, d = repmat(d,n,1); end
Z = [(d-0.16)/0.14, xi, (tw(:)-1)/3, a(:)];
end
