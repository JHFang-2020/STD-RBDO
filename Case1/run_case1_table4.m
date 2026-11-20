function run_case1_table4(seeds,nWorkers)
%RUN_CASE1_TABLE4  Reproduce Table 4 of the paper: ablation study of Case 1.
%
%   run_case1_table4()            10 independent runs per configuration
%   run_case1_table4(1:2,2)       quick check
%
%   Configurations (all with TraSUR, B = 4, budget 120 HF simulations):
%     1  ordinary Kriging                                  (SE kernel, constant mean)
%     2  Layer 2 only: SE kernel + calibrated simplified-model trend
%     3  Layer 1 only: operator constraint, constant mean
%     4  Layers 1+2, stationary discrepancy
%     5  three layers, nested scheme, inner tolerance 10 % (<= 3 batches per design)
%     6  three layers, nested scheme, inner tolerance 5 %
%     7  three layers, single loop                         (TraSUR at the current design)
%     8  three layers, single loop, design neighborhood    (proposed setting, = Table 3)
%   Results: results/main/*.mat, results/table4_case1.csv.
if nargin<1 || isempty(seeds), seeds = 1:10; end
if nargin<2 || isempty(nWorkers), nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
sp = @(a,b,ns) struct('spatial',a,'mean',b,'nonstat',ns);
cfg = {'stmsur_rb1d',sp('se','const',false),'single',0.05,Inf; ...
       'stmsur_rb1d',sp('se','lf',false),   'single',0.05,Inf; ...
       'stmsur_rb1d',sp('op','const',false),'single',0.05,Inf; ...
       'stmsur_rb1d',sp('op','lf1',false),  'single',0.05,Inf; ...
       'stmsur_rb1d',sp('op','lf1',true),   'double',0.10,3; ...
       'stmsur_rb1d',sp('op','lf1',true),   'double',0.05,Inf; ...
       'stmsur_rb1d',sp('op','lf1',true),   'single',0.05,Inf; ...
       'stmsur_rb',  sp('op','lf1',true),   'single',0.05,Inf};
jobs = struct('method',{},'B',{},'seed',{},'budget',{},'spec',{},'loop',{},'innerTol',{},'maxInner',{});
for s = seeds
    for k = 1:size(cfg,1)
        jobs(end+1) = struct('method',cfg{k,1},'B',4,'seed',s,'budget',120,'spec',cfg{k,2}, ...
            'loop',cfg{k,3},'innerTol',cfg{k,4},'maxInner',cfg{k,5}); %#ok<AGROW>
    end
end
c36_run_matrix(jobs,'main',nWorkers);
make_case1_table4(seeds);
end
