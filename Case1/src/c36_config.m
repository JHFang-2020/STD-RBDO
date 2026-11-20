function c = c36_config()
%C36_CONFIG  Case 1 (v36) problem constants and default numerical protocol.
%   G(X,xi,p,t) = R(X,p) - D(s,p),  s = (Z1(t), Z2(t), t)
%   R(X,p)      = 0.35 + X1(1-p) + X2 p            (analytical capacity)
%   D(s,p)      = A(s) - 1.8(p-c(s))^2 - 4(p-c(s))^4  (load effect, "HF snapshot")
%   X1 = d1 + m + a*zeta, X2 = d2 + m - a*zeta, a = 0.06/sqrt(2):
%   the common mode m ~ N(0,a^2) shifts G uniformly and is integrated exactly.
c = struct();
c.version     = 'case1_v36';
c.T           = 4;
c.sigmaX      = 0.06;
c.a           = 0.06/sqrt(2);
c.offset      = 0.35;
c.k           = 1.8;
c.lambda      = 4;
c.designLower = [0.4 0.4];
c.designUpper = [1.8 1.8];
c.target      = 0.01;          % allowable trajectory failure probability
% grids used by every surrogate (identical for all methods)
c.np    = 41;  c.pGrid = linspace(0,1,c.np)';      % retained spatial block of one snapshot
c.nt    = 41;  c.tGrid = linspace(0,c.T,c.nt);     % working time grid
% grids used only by the exact benchmark
c.npRef = 201; c.pRef  = linspace(0,1,c.npRef)';
c.ntRef = 401; c.tRef  = linspace(0,c.T,c.ntRef);
% state box used for space-filling initial snapshots
c.stateLower = [-3 -3 0];
c.stateUpper = [ 3  3 c.T];
end
