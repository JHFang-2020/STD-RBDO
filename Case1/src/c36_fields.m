function [A,cen] = c36_fields(S,model)
%C36_FIELDS  Demand amplitude A and moving-well centre c at load states.
%   S : M-by-3 states [Z1 Z2 t]
%   model = 'hf' : exact (high-fidelity) relations of Case 1
%         = 'lf' : simplified physics: first-order (linearised) demand
%                  amplitude about the nominal load path Z = 0, and a well
%                  wander model that neglects the Z2 cross-coupling.  It is
%                  smooth but unconservative in high-load states.
if nargin<2, model = 'hf'; end
Z1 = S(:,1); Z2 = S(:,2); t = S(:,3);
phi0 = 0.95*sin(1.1*t - pi/3);
dphi = 0.55*Z1 + 0.20*sin(pi*t/4).*Z2;
A0   = 0.64*exp(0.12*sin(pi*t/4).^2 + 0.05*t/4);
switch lower(model)
    case 'hf'
        cen = 0.5 + 0.32*tanh(phi0 + dphi);
        A   = A0.*exp(0.16*Z2 + 0.08*Z1);
    case 'lf'
        cen = 0.5 + 0.32*tanh(phi0 + 0.55*Z1);
        A   = A0.*(1 + 0.16*Z2 + 0.08*Z1);
    otherwise
        error('c36_fields:model','Unknown model %s',model);
end
end
