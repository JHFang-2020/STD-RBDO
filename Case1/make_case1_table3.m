function D = make_case1_table3(seeds)
%MAKE_CASE1_TABLE3  Table 3 of the paper from the runs in results/main.
%   For every run, the design is taken at the first iteration at which the
%   stopping rule of Eq. (22) holds at two consecutive iterations (posterior
%   s.d. of the working probability of failure <= 5 % of the target and design
%   change <= 1e-3); runs that never trigger it are taken at the budget and the
%   number of HF simulations is then reported as a lower bound "≥ mean (k/10)".
%   The probability of failure is evaluated with the exact model on the
%   independent validation population (2^17 RQMC samples).
if nargin<1 || isempty(seeds), seeds = 1:10; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
resDir = fullfile(here,'results');
T = c36_collect({'main'},'',false);                 % one row per run and iteration
save(fullfile(resDir,'T_final.mat'),'T');
T.group = T.method + "|" + T.spec;
rows = {'TraSUR (proposed)','stmsur_rb|oplf1NS'; 'SL-ABALQ','slabalq|oplf1NS'; 'EPAK','epak|oplf1NS'; ...
        'SILK','silk|oplf1NS'; 'AK-MCS-U','ukb|oplf1NS'; ...
        'Kriging + SILK','silk|seconst'; 'Kriging + AK-MCS-U','ukb|seconst'};
D = table();
fprintf('\nTable 3. Optimized designs of Case 1 (%d runs)\n',numel(seeds));
fprintf('%-34s %-18s %-20s %-20s %s\n','Method','HF simulations','Design (d1, d2)','Objective','Pf (validation)');
R = load(fullfile(resDir,'reference_rows.mat')); R = R.R;
for k = 1:height(R)
    if R.pf(k) >= 0.1, pfs = sprintf('%.1f%%',100*R.pf(k)); elseif R.pf(k) >= 0.02, pfs = sprintf('%.2f%%',100*R.pf(k));
    else, pfs = sprintf('%.3f%%',100*R.pf(k)); end
    fprintf('%-34s %-18s (%.4f, %.4f)       %-20s %s\n',R.name(k), ...
        R.hf(k),R.d1(k),R.d2(k),sprintf('%.4f',R.C(k)),pfs);
end
for k = 1:size(rows,1)
    X = nan(numel(seeds),5); trig = false(numel(seeds),1);
    for i = 1:numel(seeds)
        Rk = sortrows(T(T.group==rows{k,2} & T.seed==seeds(i) & T.B==4,:),'nHF');
        if isempty(Rk), continue; end
        dmv = [Inf; max(abs(diff([Rk.d1 Rk.d2],1,1)),[],2)];
        okS = Rk.relSd <= 0.05 & dmv <= 1e-3;
        js = find(okS(2:end) & okS(1:end-1),1) + 1;
        trig(i) = ~isempty(js); if ~trig(i), js = height(Rk); end
        X(i,:) = [Rk.nHF(js) Rk.d1(js) Rk.d2(js) Rk.d1(js)+Rk.d2(js) Rk.pfTrue(js)];
    end
    ok = all(isfinite(X),2); X = X(ok,:); trig = trig(ok); n = size(X,1);
    if n == 0, continue; end
    if all(trig), hf = sprintf('%.1f ± %.1f',mean(X(:,1)),std(X(:,1)));
    else, hf = sprintf('≥%.1f (%d/%d)',mean(X(:,1)),nnz(trig),n); end
    fprintf('%-34s %-18s (%.4f, %.4f)       %.4f ± %.4f      %.3f%% ± %.3f%%\n',rows{k,1},hf, ...
        mean(X(:,2)),mean(X(:,3)),mean(X(:,4)),std(X(:,4)),100*mean(X(:,5)),100*std(X(:,5)));
    D = [D; table(string(rows{k,1}),n,nnz(trig),mean(X(:,1)),std(X(:,1)),mean(X(:,2)),mean(X(:,3)), ...
        mean(X(:,4)),std(X(:,4)),mean(X(:,5)),std(X(:,5)),min(X(:,5)),max(X(:,5)), ...
        'VariableNames',{'method','nRuns','nStopped','HFmean','HFstd','d1','d2','Cmean','Cstd', ...
        'pfValMean','pfValStd','pfValMin','pfValMax'})]; %#ok<AGROW>
end
writetable(D,fullfile(resDir,'table3_case1.csv'));
fprintf('written: %s\n',fullfile(resDir,'table3_case1.csv'));
end
