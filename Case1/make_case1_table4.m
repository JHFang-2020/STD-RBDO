function D = make_case1_table4(seeds)
%MAKE_CASE1_TABLE4  Table 4 of the paper (ablation) from the runs in results/main.
%   N_M / N_H: HF simulations after which medium accuracy (eps_P <= 5 % and
%   eps_C <= 0.3 %) / high accuracy (eps_P <= 2 % and eps_C <= 0.1 %) hold at
%   every later iteration, with eps_P and eps_C measured against the exact
%   optimum on the run's own working population (Eq. (23)); runs that never get
%   there are counted as 120 and the mean is shown as "≥ mean (k/10)".
%   eps_P(n): mean |eps_P| at the first iteration with at least n HF simulations.
if nargin<1 || isempty(seeds), seeds = 1:10; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
resDir = fullfile(here,'results');
T = c36_collect({'main'},'',false);                 % incremental: unchanged runs are reused
save(fullfile(resDir,'T_final.mat'),'T');
T.group = T.method + "|" + T.spec + "|" + T.loop;
rows = {'1','stmsur_rb1d|seconst|single'; '2','stmsur_rb1d|self|single'; '3','stmsur_rb1d|opconst|single'; ...
        '4','stmsur_rb1d|oplf1|single'; '5','stmsur_rb1d|oplf1NS|double_m3'; '6','stmsur_rb1d|oplf1NS|double'; ...
        '7','stmsur_rb1d|oplf1NS|single'; '8 (Proposed)','stmsur_rb|oplf1NS|single'};
tol = [0.05 0.003; 0.02 0.001];
D = table();
fprintf('\nTable 4. Ablation study of Case 1 (means of %d runs)\n',numel(seeds));
fprintf('%-14s %-16s %-16s %-14s %s\n','Config.','N_M','N_H','eps_P(n=36)','eps_P(n=120)');
for k = 1:size(rows,1)
    N = 120*ones(numel(seeds),2); E = nan(numel(seeds),2); ok = false(numel(seeds),1);
    for i = 1:numel(seeds)
        R = sortrows(T(T.group==rows{k,2} & T.seed==seeds(i) & T.B==4,:),'nHF');
        if isempty(R), continue; end
        ok(i) = true; reached = false(1,2);
        for q = 1:2
            g = abs(R.feasErrW) <= tol(q,1) & abs(R.costGapW) <= tol(q,2);
            last = find(~g,1,'last');
            if isempty(last), N(i,q) = R.nHF(1); reached(q) = true;
            elseif last < height(R), N(i,q) = R.nHF(last+1); reached(q) = true; end
        end
        N(i,3:4) = reached;
        E(i,1) = abs(R.feasErrW(find(R.nHF >= 36,1))); E(i,2) = abs(R.feasErrW(find(R.nHF >= 120,1)));
    end
    N = N(ok,:); E = E(ok,:); n = size(N,1);
    if n == 0, continue; end
    s = cell(1,2);
    for q = 1:2
        kq = nnz(N(:,2+q));
        if kq == n, s{q} = sprintf('%.1f',mean(N(:,q))); else, s{q} = sprintf('≥%.1f (%d/%d)',mean(N(:,q)),kq,n); end
    end
    fprintf('%-14s %-16s %-16s %-14s %.2f%%\n',rows{k,1},s{1},s{2},sprintf('%.1f%%',100*mean(E(:,1))),100*mean(E(:,2)));
    D = [D; table(string(rows{k,1}),n,mean(N(:,1)),nnz(N(:,3)),mean(N(:,2)),nnz(N(:,4)),mean(E(:,1)),mean(E(:,2)), ...
        'VariableNames',{'config','nRuns','NM','NMreached','NH','NHreached','epsP36','epsP120'})]; %#ok<AGROW>
end
writetable(D,fullfile(resDir,'table4_case1.csv'));
fprintf('written: %s\n',fullfile(resDir,'table4_case1.csv'));
end
