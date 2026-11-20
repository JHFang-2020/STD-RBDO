function K = c36_ckern(cm,S1,S2)
%C36_CKERN  Prior covariance of one discrepancy component.
%   K = sf2 * w(S1) w(S2)' .* k_M52(S1,S2)
%   w(s)^2 = 1 + gamma*(Z1^2 + Z2^2): physics-informed non-stationary
%   amplitude (the linearised-load model error grows with the load
%   perturbation |Z|); gamma = 0 recovers the stationary kernel.
K = cm.sf2*c36_kernel(S1,S2,cm.ell,'m52');
if isfield(cm,'gamma') && cm.gamma > 0
    K = K.*(c36_wamp(S1,cm.gamma)*c36_wamp(S2,cm.gamma)');
end
end
