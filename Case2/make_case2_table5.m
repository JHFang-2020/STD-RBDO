function D = make_case2_table5(seeds,nWorkers)
%MAKE_CASE2_TABLE5  Table 5 of the paper from the runs in results/main.
%   For every run, the design is taken at the first iteration at which the
%   stopping rule of Eq. (22) holds at two consecutive iterations (posterior
%   s.d. of the working probability of failure <= 5 % of the target and
%   max_i |d_i,k - d_i,k-1| <= 1e-3 m); runs that never trigger it are taken
%   at the budget (80 HF simulations) and the number of HF simulations is
%   then reported as a lower bound "≥ mean (k/n)".
%   The probability of failure of every stopped design is evaluated on the
%   independent validation population (16,384 histories, full-grid HF,
%   b36_validate). Validations are cached in results/C2_val_cache.mat; only
%   designs that are not in the cache are validated (nWorkers parallel workers).
%   The first three rows are read from results/reference_rows.mat
%   (case2_reference_rows, computed first if missing).
%   Writes results/table5_case2.csv (table) and results/table5_case2_runs.csv
%   (stopped design of every run).
if nargin<1 || isempty(seeds), seeds = 1:10; end
if nargin<2 || isempty(nWorkers), nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
resDir = fullfile(here,'results');
T = b36_collect({'main'});                           % one row per run and iteration
save(fullfile(resDir,'T_final.mat'),'T');
T.group = T.method + "|" + T.spec + "|" + T.loop + "|B" + string(T.B);
rows = {'TraSUR (proposed)','stmsur|EfieldNS|single|B4'; 'SL-ABALQ','slabalq|EfieldNS|single|B4'; ...
        'EPAK','epak|EfieldNS|single|B4'; 'SILK','silk|EfieldNS|single|B4'; 'AK-MCS-U','ukb|EfieldNS|single|B4'; ...
        'Layer 1 only (ablation)','stmsur|qs|single|B4'};
% ---- stopped design of every run
L = table();
for k = 1:size(rows,1)
    for i = 1:numel(seeds)
        R = sortrows(T(T.group==rows{k,2} & T.seed==seeds(i),:),'nHF');
        if isempty(R), continue; end
        dmv = [Inf; max(abs(diff([R.d1 R.d2 R.d3],1,1)),[],2)];
        okS = R.relSd <= 0.05 & dmv <= 1e-3;
        js = find(okS(2:end) & okS(1:end-1),1) + 1;
        trig = ~isempty(js); if ~trig, js = height(R); end
        L = [L; table(string(rows{k,1}),string(rows{k,2}),seeds(i),trig,R.nHF(js),[R.d1(js) R.d2(js) R.d3(js)], ...
            'VariableNames',{'method','group','seed','stopped','nHF','d'})]; %#ok<AGROW>
    end
end
if isempty(L), error('make_case2_table5:noRuns','no runs found in %s',fullfile(resDir,'main')); end
% ---- validation of the stopped designs (cached)
cf = fullfile(resDir,'C2_val_cache.mat');
if isfile(cf), s = load(cf); Vc = s.Vc; else, Vc = emptyCache(); end
need = unique(L.d,'rows','stable');
have = ~isnan(lookupPf(Vc,need));
if any(~have)
    todo = need(~have,:);
    fprintf('validating %d designs on the validation population\n',size(todo,1));
    pool = b36_pool(nWorkers); nb = 16; if ~isempty(pool), nb = 2*pool.NumWorkers; end
    for b0 = 1:nb:size(todo,1)
        Vn = b36_validate(todo(b0:min(end,b0+nb-1),:));
        Vc = [Vc; Vn(:,Vc.Properties.VariableNames)]; %#ok<AGROW>
        save(cf,'Vc');
    end
end
L.pfVal = lookupPf(Vc,L.d);
Lw = L; Lw.d1 = L.d(:,1); Lw.d2 = L.d(:,2); Lw.d3 = L.d(:,3); Lw.d = [];
writetable(Lw,fullfile(resDir,'table5_case2_runs.csv'));
% ---- table
fr = fullfile(resDir,'reference_rows.mat');
if ~isfile(fr), case2_reference_rows(nWorkers); end
Rf = load(fr); Rf = Rf.R;
fprintf('\nTable 5. Optimized designs at stopping in Case 2 (%d runs)\n',numel(seeds));
fprintf('%-34s %-16s %-28s %-20s %s\n','Method','HF simulations','Design (d1, d2, d3)','Objective','Pf (validation)');
D = table();
for k = 1:height(Rf)
    if Rf.pf(k) >= 0.1, pfs = sprintf('%.1f%%',100*Rf.pf(k)); else, pfs = sprintf('%.3f%%',100*Rf.pf(k)); end
    fprintf('%-34s %-16s (%.4f, %.4f, %.5f)     %-20s %s\n',Rf.name(k),Rf.hf(k),Rf.d1(k),Rf.d2(k),Rf.d3(k), ...
        sprintf('%.5f',Rf.C(k)),pfs);
    D = [D; table(Rf.name(k),NaN,NaN,NaN,NaN,Rf.d1(k),Rf.d2(k),Rf.d3(k),Rf.C(k),NaN,Rf.pf(k),NaN,NaN,NaN, ...
        'VariableNames',varNames())]; %#ok<AGROW>
end
for k = 1:size(rows,1)
    Lk = L(L.method==rows{k,1},:); n = height(Lk);
    if n == 0, continue; end
    X = [Lk.nHF Lk.d sum(Lk.d,2) Lk.pfVal]; trig = Lk.stopped;
    if all(trig), hf = sprintf('%.1f ± %.1f',mean(X(:,1)),std(X(:,1)));
    else, hf = sprintf('≥%.1f (%d/%d)',mean(X(:,1)),nnz(trig),n); end
    fprintf('%-34s %-16s (%.4f, %.4f, %.5f)     %.5f ± %.5f    %.3f%% ± %.3f%%\n',rows{k,1},hf, ...
        mean(X(:,2)),mean(X(:,3)),mean(X(:,4)),mean(X(:,5)),std(X(:,5)),100*mean(X(:,6)),100*std(X(:,6)));
    D = [D; table(string(rows{k,1}),n,nnz(trig),mean(X(:,1)),std(X(:,1)),mean(X(:,2)),mean(X(:,3)),mean(X(:,4)), ...
        mean(X(:,5)),std(X(:,5)),mean(X(:,6)),std(X(:,6)),min(X(:,6)),max(X(:,6)),'VariableNames',varNames())]; %#ok<AGROW>
end
writetable(D,fullfile(resDir,'table5_case2.csv'));
fprintf('written: %s\n',fullfile(resDir,'table5_case2.csv'));
end

function v = varNames()
v = {'method','nRuns','nStopped','HFmean','HFstd','d1','d2','d3','Cmean','Cstd','pfValMean','pfValStd','pfValMin','pfValMax'};
end

function p = lookupPf(Vc,dd)
[tf,loc] = ismembertol(dd,[Vc.d1 Vc.d2 Vc.d3],1e-12,'ByRows',true,'DataScale',1);
p = nan(size(dd,1),1); p(tf) = Vc.pfVal(loc(tf));
end

function Vc = emptyCache()
z = zeros(0,1);
Vc = table(z,z,z,z,z,z,z,z,false(0,1),z,'VariableNames', ...
    {'d1','d2','d3','cost','pfVal','pfValSE','nFailVal','qVal','certified','nHF'});
end
