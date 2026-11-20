function out = b36_ddo()
%B36_DDO  Deterministic design optimum of Case 2 and its true failure probability.
%   Nominal inputs (xi = 0: nominal geometry, modulus and load histories):
%   min d1 + d2 + d3  s.t.  max_{x,t}|sigma(d,0)| <= 160 MPa  (HF, full grid).
%   The failure probability of the resulting design is evaluated on the
%   independent validation bank (b36_pf_val).
c = b36_config(); xi0 = zeros(1,17);
pk = @(d) getfield(b36_hf(d(:)',xi0,'full'),'peak');
% the stress constraint is decreasing in every d_i: start from the
% structure of the probabilistic optimum and let SQP move all three heights
opt = optimoptions('fmincon','Algorithm','sqp','Display','iter','FiniteDifferenceStepSize',1e-4, ...
    'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-6,'StepTolerance',1e-7);
d0 = [0.16 0.30 0.17];
d = fmincon(@(d) sum(d),d0,[],[],[],[],c.lower,c.upper,@(d) deal(pk(d)/c.allow - 1,[]),opt);
V = b36_pf_val(d);
out = struct('d',d,'cost',sum(d),'peakNominal',pk(d),'pf',V.pfVal,'pfSE',V.pfValSE,'nFail',V.nFailVal,'certified',V.certified,'nHF',V.nHF);
fprintf('DDO: d = (%.5f, %.5f, %.5f), C = %.5f, nominal peak %.2f MPa, Pf(validation) = %.3f%% (%d failures, certified %d)\n', ...
    d,sum(d),out.peakNominal,100*out.pf,out.nFail,out.certified);
save(fullfile(fileparts(mfilename('fullpath')),'..','results','C2_ddo.mat'),'out');
end
