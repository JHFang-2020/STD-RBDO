function run_case2_table5(seeds,nWorkers)
%RUN_CASE2_TABLE5  Reproduce Table 5 of the paper: optimized designs of Case 2.
%
%   run_case2_table5()            10 independent runs per configuration (as in the paper)
%   run_case2_table5(1:2,2)       quick check: runs 1-2 only, 2 parallel workers
%
%   Six configurations are run with batch size B = 4 and a budget of 80 HF
%   simulations (the 16 initial ones included):
%     TraSUR (proposed, with the design neighbourhood), SL-ABALQ, EPAK, SILK
%     and AK-MCS-U on the three-layer PI-Kriging surrogate, and the
%     "Layer 1 only" ablation (TraSUR with the exact quasi-static moment as
%     the only physics and a stationary discrepancy).
%   Every run continues to the budget; the common stopping rule of Eq. (22)
%   is applied afterwards to the recorded iterations (make_case2_table5).
%   Finished runs are skipped, so the script can be restarted.
%   Results: results/main/*.mat (one file per run), results/table5_case2.csv.
if nargin<1 || isempty(seeds), seeds = 1:10; end
if nargin<2 || isempty(nWorkers), nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
piK = struct('mean','lf','nonstat',true,'lf','Efield');    % Layers 1-3: exact statics + two-mode simplified model + non-stationary discrepancy
qsK = struct('mean','lf','nonstat',false,'lf','qs');       % Layer 1 only: exact quasi-static moment, stationary discrepancy
cfg = {'stmsur',piK; 'slabalq',piK; 'epak',piK; 'silk',piK; 'ukb',piK; 'stmsur',qsK};
jobs = struct('method',{},'B',{},'seed',{},'budget',{},'spec',{},'loop',{},'innerTol',{},'maxInner',{});
for s = seeds
    for k = 1:size(cfg,1)
        jobs(end+1) = struct('method',cfg{k,1},'B',4,'seed',s,'budget',80,'spec',cfg{k,2}, ...
            'loop','single','innerTol',0.05,'maxInner',Inf); %#ok<AGROW>
    end
end
b36_run_matrix(jobs,'main',nWorkers);
case2_reference_rows(nWorkers);
make_case2_table5(seeds,nWorkers);
end
