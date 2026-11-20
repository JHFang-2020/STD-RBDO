function run_case1_table3(seeds,nWorkers)
%RUN_CASE1_TABLE3  Reproduce Table 3 of the paper: optimized designs of Case 1.
%
%   run_case1_table3()            10 independent runs per method (as in the paper)
%   run_case1_table3(1:2,2)       quick check: runs 1-2 only, 2 parallel workers
%
%   Seven configurations are run with batch size B = 4 and a budget of 120 HF
%   simulations (the 12 initial ones included):
%     TraSUR (proposed), SL-ABALQ, EPAK, SILK and AK-MCS-U on the three-layer
%     PI-Kriging surrogate, and SILK and AK-MCS-U on ordinary Kriging.
%   Every run continues to the budget; the common stopping rule of Eq. (22) is
%   applied afterwards to the recorded iterations (make_case1_table3).
%   Finished runs are skipped, so the script can be restarted.
%   Results: results/main/*.mat (one file per run), results/table3_case1.csv.
if nargin<1 || isempty(seeds), seeds = 1:10; end
if nargin<2 || isempty(nWorkers), nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
piK = struct('spatial','op','mean','lf1','nonstat',true);   % three-layer PI-Kriging (Layers 1-3)
okK = struct('spatial','se','mean','const');                % ordinary Kriging
cfg = {'stmsur_rb',piK; 'slabalq',piK; 'epak',piK; 'silk',piK; 'ukb',piK; 'silk',okK; 'ukb',okK};
jobs = struct('method',{},'B',{},'seed',{},'budget',{},'spec',{},'loop',{},'innerTol',{},'maxInner',{});
for s = seeds
    for k = 1:size(cfg,1)
        jobs(end+1) = struct('method',cfg{k,1},'B',4,'seed',s,'budget',120,'spec',cfg{k,2}, ...
            'loop','single','innerTol',0.05,'maxInner',Inf); %#ok<AGROW>
    end
end
c36_run_matrix(jobs,'main',nWorkers);
case1_reference_rows();
make_case1_table3(seeds);
end
