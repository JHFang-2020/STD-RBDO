function c36_run_matrix(jobs,outName,nWorkers)
%C36_RUN_MATRIX  Run a list of job configs in parallel; skip finished ones.
%   jobs : struct array with fields method,B,seed,budget,spec,(loop)
if nargin<3, nWorkers = 9; end
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here,'..','results',outName);
if ~isfolder(outDir), mkdir(outDir); end
seeds = unique([jobs.seed]);
c36_truth('validation');
for s = seeds, c36_truth('work',s); end
% skip completed jobs
todo = true(numel(jobs),1);
for k = 1:numel(jobs)
    j = jobs(k); j.resDir = outDir; j.verbose = false;
    lab = labelOf(j);
    if isfile(fullfile(outDir,[lab '.mat'])), todo(k) = false; end
end
jobs = jobs(todo);
fprintf('%d jobs to run in %s\n',numel(jobs),outDir);
if isempty(jobs), return; end
% longest jobs first
cost = arrayfun(@(j) jobCost(j),jobs); [~,o] = sort(cost,'descend'); jobs = jobs(o);
c36_pool(nWorkers);
logf = fullfile(outDir,'run_log.txt');
parfor k = 1:numel(jobs)
    j = jobs(k); j.resDir = outDir; j.verbose = false;
    t0 = tic;
    try
        r = c36_run(j);
        msg = sprintf('%s done %.0fs final nHF=%d C=%.4f Pf=%.4f%%',r.cfg.label,toc(t0),r.final.nHF,r.final.cost,100*r.final.pfTrue);
    catch err
        msg = sprintf('%s FAILED: %s',labelOf(j),err.message);
    end
    fid = fopen(logf,'a'); fprintf(fid,'%s\n',msg); fclose(fid);
end
end

function c = jobCost(j)
w = struct('stmsur_mc',4,'pointsur',4,'stmsur',2,'stmsur_1d',2,'rvsur',2.5,'rvsur_1d',2,'tmp',1.5);
if isfield(w,j.method), c = w.(j.method); else, c = 1; end
c = c*j.budget/max(j.B,1)^0.7;
end

function lab = labelOf(j)
ns = ''; if isfield(j.spec,'nonstat') && j.spec.nonstat, ns = 'NS'; end
lp = ''; if isfield(j,'loop') && strcmp(j.loop,'double'), lp = '_DL'; end
if ~isempty(lp) && isfield(j,'maxInner') && isfinite(j.maxInner), lp = sprintf('_DLm%d',j.maxInner); end
lab = sprintf('%s_%s%s%s_B%d_s%03d%s',j.method,j.spec.spatial,j.spec.mean,ns,j.B,j.seed,lp);
end
