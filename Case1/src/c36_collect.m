function T = c36_collect(dirNames,outCsv,skipVal)
%C36_COLLECT  Long table (one row per run and checkpoint) from saved runs.
if ischar(dirNames), dirNames = {dirNames}; end
here = fileparts(mfilename('fullpath'));
c = c36_config();
rows = {}; wcache = containers.Map('KeyType','double','ValueType','any');
% incremental: reuse rows of a previously collected table for unchanged runs
cacheFile = fullfile(here,'..','results','T_final.mat'); Tprev = []; tPrev = 0;
if nargin>=3 && skipVal, cacheFile = ''; end
if ~isempty(cacheFile) && isfile(cacheFile)
    sp = load(cacheFile,'T'); Tprev = sp.T; di = dir(cacheFile); tPrev = di.datenum;
    if ismember('group',Tprev.Properties.VariableNames), Tprev.group = []; end
end
for q = 1:numel(dirNames)
    d = dir(fullfile(here,'..','results',dirNames{q},'*.mat'));
    for k = 1:numel(d)
        lab = erase(d(k).name,'.mat');
        if ~isempty(Tprev) && d(k).datenum < tPrev
            Rk = Tprev(Tprev.label==lab & Tprev.dir==dirNames{q},:);
            if ~isempty(Rk), rows{end+1} = Rk; continue; end %#ok<AGROW>
        end
        s = load(fullfile(d(k).folder,d(k).name)); r = s.res; h = r.hist; cfg = r.cfg;
        % population-attributable reference: exact optimum on the run's own
        % working population (same histories, same grids) -> isolates the
        % surrogate/acquisition error from the Monte Carlo error of the population
        if ~isKey(wcache,cfg.seed)
            W = c36_truth('work',cfg.seed);
            solW = c36_solve_rbdo(W.Etrue,c.pGrid,W.bank.zeta,c.target);
            wcache(cfg.seed) = struct('E',W.Etrue,'zeta',W.bank.zeta,'CW',solW.cost,'dW',solW.d);
        end
        wc = wcache(cfg.seed);
        dd = reshape([h.d],2,[])';
        % validation-bank scoring (post hoc; identical function to the in-run version)
        if nargin<3, skipVal = false; end
        if any(isnan([h.pfTrue])) && ~skipVal
            if ~exist('VALc','var'), VALc = c36_truth('validation'); VALc.Ed = double(VALc.E); end
            for i = 1:numel(h)
                [h(i).pfTrue,~,h(i).pfTrueSE] = c36_pf(VALc.Ed,VALc.pGrid,dd(i,:),VALc.zeta);
            end
        end
        pfW = arrayfun(@(i) c36_pf(wc.E,c.pGrid,dd(i,:),wc.zeta),(1:size(dd,1))');
        ns = isfield(cfg.spec,'nonstat') && cfg.spec.nonstat;
        loop = 'single'; if isfield(cfg,'loop'), loop = cfg.loop; end
        if strcmp(loop,'double') && isfield(cfg,'maxInner') && isfinite(cfg.maxInner)
            loop = sprintf('double_m%d',cfg.maxInner);
        end
        n = numel(h);
        tab = table(repmat(string(dirNames{q}),n,1),repmat(string(cfg.label),n,1),repmat(string(cfg.method),n,1), ...
            repmat(string([cfg.spec.spatial cfg.spec.mean ternary(ns,'NS','')]),n,1),repmat(string(loop),n,1), ...
            repmat(cfg.B,n,1),repmat(cfg.seed,n,1),(1:n)',[h.nHF]',[h.cost]',[h.costGap]',[h.pfTrue]',[h.pfTrueSE]', ...
            [h.pfSur]',[h.pfRelErrStar]',[h.hRMSE]',[h.hRMSEcrit]',[h.misclass]',[h.epsHat]',[h.relSd]', ...
            [h.elapsed]',[h.tFit]',[h.tPop]',[h.tOpt]',[h.tEv]',[h.tAcq]', ...
            reshape([h.d],2,[])', ...
            'VariableNames',{'dir','label','method','spec','loop','B','seed','iter','nHF','cost','costGap','pfTrue','pfTrueSE', ...
            'pfSur','pfRelErrStar','hRMSE','hRMSEcrit','misclass','epsHat','relSd','elapsed','tFit','tPop','tOpt','tEv','tAcq','d'});
        tab.pfTrueW = pfW;
        tab.feasErrW = (pfW - c.target)/c.target;
        tab.costGapW = ([h.cost]' - wc.CW)/wc.CW;
        tab.CW = repmat(wc.CW,n,1);
        % calibration of the posterior s.d. of the working risk at the incumbent
        tab.sdPf = [h.sdPf]';
        tab.zRisk = (tab.pfSur - pfW)./max(tab.sdPf,1e-12);
        tab.d1 = tab.d(:,1); tab.d2 = tab.d(:,2); tab.d = [];
        tab.feasErr = (tab.pfTrue - 0.01)/0.01;
        if ~isempty(Tprev), tab = tab(:,Tprev.Properties.VariableNames); end
        rows{end+1} = tab; %#ok<AGROW>
    end
end
T = vertcat(rows{:});
if nargin>1 && ~isempty(outCsv), writetable(T,outCsv); end
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
