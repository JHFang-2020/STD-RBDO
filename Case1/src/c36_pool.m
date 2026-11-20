function pool = c36_pool(n)
%C36_POOL  Process pool with n single-threaded workers.
%   Without the Parallel Computing Toolbox, no pool is started and the parfor
%   loops of c36_run_matrix run serially.
pool = [];
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox'), return; end
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers == n, return; end
if ~isempty(pool), delete(pool); end
cl = parcluster('Processes');
cl.NumWorkers = max(cl.NumWorkers,n);
pool = parpool(cl,n);
end
