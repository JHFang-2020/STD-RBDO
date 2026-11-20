function T = b36_collect(dirNames,outCsv)
%B36_COLLECT  Long table (one row per run and checkpoint) of Case 2 runs.
%   Errors are w.r.t. the exact optimum of the SAME working population
%   (b36_truth_ref): surrogate/acquisition-attributable.
if ischar(dirNames), dirNames = {dirNames}; end
here = fileparts(mfilename('fullpath'));
rows = {};
for q = 1:numel(dirNames)
    dd = dir(fullfile(here,'..','results',dirNames{q},'*.mat'));
    for k = 1:numel(dd)
        s = load(fullfile(dd(k).folder,dd(k).name)); r = s.res; h = r.hist; cfg = r.cfg;
        n = numel(h);
        sp = cfg.spec.lf; if cfg.spec.nonstat, sp = [sp 'NS']; end
        loop = cfg.loop; if strcmp(loop,'double') && isfinite(cfg.maxInner), loop = sprintf('double_m%d',cfg.maxInner); end
        if strcmp(cfg.loop,'double') && ~isfinite(cfg.maxInner), loop = sprintf('double_tol%02d',round(100*cfg.innerTol)); end
        D = reshape([h.d],3,[])';
        tab = table(repmat(string(dirNames{q}),n,1),repmat(string(cfg.label),n,1),repmat(string(cfg.method),n,1), ...
            repmat(string(sp),n,1),repmat(string(loop),n,1),repmat(cfg.B,n,1),repmat(cfg.seed,n,1),(1:n)', ...
            [h.nHF]',[h.cost]',[h.costGapW]',[h.qTrue]',[h.qErr]',[h.nFailW]',[h.pfTrueW]',[h.feasErrW]', ...
            [h.pfSur]',[h.pfPlug]',[h.sdPf]',[h.relSd]',[h.pfSurStar]',[h.pfRelErrStar]',[h.hRMSEcrit]', ...
            [h.misclassStar]',[h.nAmbiguous]',[h.nActive]',[h.truthExact]',[h.elapsed]',[h.tFit]',[h.tOpt]',[h.tEv]',[h.tAcq]', ...
            D(:,1),D(:,2),D(:,3), ...
            'VariableNames',{'dir','label','method','spec','loop','B','seed','iter','nHF','cost','costGap','qTrue','qErr', ...
            'nFailW','pfTrueW','feasErrW','pfSur','pfPlug','sdPf','relSd','pfSurStar','pfRelErrStar','hRMSEcrit', ...
            'misclassStar','nAmbiguous','nActive','truthExact','elapsed','tFit','tOpt','tEv','tAcq','d1','d2','d3'});
        tab.zRisk = (tab.pfSur - tab.pfTrueW)./max(tab.sdPf,1e-12);
        rows{end+1} = tab; %#ok<AGROW>
    end
end
T = vertcat(rows{:});
if nargin>1 && ~isempty(outCsv), writetable(T,outCsv); end
end
