function b36_run_matrix(jobs,outName,nWorkers)
%B36_RUN_MATRIX  Run a list of Case 2 jobs in parallel (one run per worker);
%   finished runs are skipped.
if nargin<3, nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here,'..','results',outName);
if ~isfolder(outDir), mkdir(outDir); end
if ~isfile(fullfile(here,'..','results','truth','truth_ref.mat')), b36_pool(nWorkers); end   % built on the same pool
b36_truth_ref();                                   % must exist before the workers start
todo = true(numel(jobs),1); labs = cell(numel(jobs),1);
for k = 1:numel(jobs)
    labs{k} = labelOf(jobs(k));
    if isfile(fullfile(outDir,[labs{k} '.mat'])), todo(k) = false; end
end
jobs = jobs(todo); labs = labs(todo);
[~,u] = unique(string(labs),'stable'); jobs = jobs(u); labs = labs(u);   % a configuration shared by two parts runs once
fprintf('%d jobs to run in %s\n',numel(jobs),outDir);
if isempty(jobs), return; end
% seed by seed (interim analyses possible), longest jobs first within a seed
cost = arrayfun(@(j) (1 + strcmp(j.method,'stmsur'))*j.budget/max(j.B,1)^0.7,jobs);
[~,o] = sortrows([[jobs.seed]' -cost(:)]); jobs = jobs(o); labs = labs(o);
b36_pool(nWorkers);
logf = fullfile(outDir,'run_log.txt');
parfor k = 1:numel(jobs)
    j = jobs(k); j.resDir = outDir; j.verbose = false; j.label = labs{k};
    t0 = tic;
    try
        r = b36_run(j);
        msg = sprintf('%s done %.0fs final nHF=%d C=%.5f gap=%+.4f%% qErr=%+.3f%%',labs{k},toc(t0),r.final.nHF,r.final.cost,100*r.final.costGapW,100*r.final.qErr);
    catch err
        msg = sprintf('%s FAILED: %s (%s line %d)',labs{k},err.message,err.stack(1).name,err.stack(1).line);
    end
    fid = fopen(logf,'a'); fprintf(fid,'%s %s\n',char(datetime('now','Format','HH:mm:ss')),msg); fclose(fid);
end
end

function lab = labelOf(j)
sp = j.spec.lf; if j.spec.nonstat, sp = [sp 'NS']; end
lp = ''; if strcmp(j.loop,'double'), lp = '_DL'; if isfinite(j.maxInner), lp = sprintf('_DLm%d',j.maxInner); end, end
lab = sprintf('%s_%s_B%d_s%03d%s',j.method,sp,j.B,j.seed,lp);
end
