function res = c36_run(cfg)
%C36_RUN  One single-loop STD-RBDO run (one method, batch size, seed).
%   Loop:  fit surrogate -> solve the frozen surrogate RBDO globally ->
%          record checkpoint -> select a batch of snapshot states ->
%          evaluate them (HF) -> append.  Identical for every method; only
%          the acquisition (cfg.method) and the surrogate spec differ.
%   cfg fields: method, B, seed, budget, N0, spec, resDir, label
c = c36_config();
cfg = defaults(cfg);
VAL = c36_truth('validation',[],cfg.truthDir);
WRK = c36_truth('work',cfg.seed,cfg.truthDir);
bank = WRK.bank;
dstar = VAL.dstar; Cstar = VAL.Cstar;
if cfg.evalVal, Eval = double(VAL.E); valP = VAL.pGrid; valZ = VAL.zeta; end
clear VAL                       % validation scoring is done post hoc (c36_collect) by default
% ---------- initial snapshots (identical for all methods of a seed)
S = c36_lhs_states(cfg.N0,1000+cfg.seed);
Y = c36_loadeffect(S,c.pGrid,'hf');
nHF = size(S,1);
model = []; hist = struct([]); it = 0; tStart = tic;
hTrueStar = c36_margin(WRK.Etrue,c.pGrid,dstar,bank.zeta);
pfTrueWstar = mean(0.5*erfc(hTrueStar/(sqrt(2)*c.a)));
mN = c.a*bank.um;
dPrev = []; dCur = []; innerDone = false; nOuter = 0; nInner = 0;
while true
    it = it + 1;
    t1 = tic; model = c36_cgp_fit(S,Y,cfg.spec,model); tFit = toc(t1);
    t1 = tic; pop = c36_popeval(model,bank); tPop = toc(t1);
    t1 = tic;
    sol = c36_solve_rbdo(pop.E,c.pGrid,bank.zeta,c.target);
    if ~sol.feasible
        sol.d = c.designUpper; sol.cost = sum(sol.d); sol.h = 0; sol.m0 = 1.8;
        sol.pf = c36_pf(pop.E,c.pGrid,sol.d,bank.zeta);
    end
    dk = sol.d; tOpt = toc(t1);
    if strcmp(cfg.loop,'double')
        % nested coupling with data reuse: the design is only moved to the
        % surrogate optimum once the inner reliability analysis at the
        % current design has converged (relSd <= innerTol)
        if isempty(dCur) || innerDone || nInner >= cfg.maxInner
            dCur = dk; innerDone = false; nOuter = nOuter + 1; nInner = 0;
        end
        nInner = nInner + 1;
        dk = dCur;
        sol.d = dk; sol.h = (dk(1)-dk(2))/2; sol.m0 = mean(dk);
    end
    % ---------- events at the incumbent (shared diagnostics for all methods)
    t1 = tic;
    P0 = c36_pointfield(pop,bank,dk);
    ev0 = c36_events(model,bank,P0,struct('K',cfg.K,'Emax',cfg.Emax));
    F0 = c36_rvfun(model,bank,P0,struct('K',cfg.K,'Tmax',cfg.Tmax));
    tEv = toc(t1);
    % ---------- checkpoint metrics (truth never fed back)
    rec = struct();
    rec.iter = it; rec.nHF = nHF; rec.d = dk; rec.cost = sum(dk);
    if cfg.evalVal, [rec.pfTrue,~,rec.pfTrueSE] = c36_pf(Eval,valP,dk,valZ);
    else, rec.pfTrue = NaN; rec.pfTrueSE = NaN; end
    rec.pfSur = c36_pf(pop.E,c.pGrid,dk,bank.zeta);
    hS = c36_margin(pop.E,c.pGrid,dstar,bank.zeta);
    rec.pfSurStar = mean(0.5*erfc(hS/(sqrt(2)*c.a)));
    rec.pfTrueWstar = pfTrueWstar;
    rec.pfRelErrStar = (rec.pfSurStar - pfTrueWstar)/pfTrueWstar;
    rec.hRMSE = rms(hS - hTrueStar);
    crit = abs(hTrueStar) < 3*c.a;
    rec.hRMSEcrit = rms(hS(crit) - hTrueStar(crit));
    rec.misclass = sum(sign(hS + mN) ~= sign(hTrueStar + mN));
    pin = ev0.pi;
    rec.expMisclass = sum(min(pin,1-pin));
    rec.nPredFail = sum(pin > 0.5);
    rec.epsHat = rec.expMisclass/max(rec.nPredFail,1);
    rec.nAmbiguous = numel(ev0.amb);
    rec.sdPf = F0.sdPf;                       % posterior s.d. of the working risk
    rec.relSd = F0.sdPf/c.target;
    rec.pfSurAtDk = c36_pf(pop.E,c.pGrid,dk,bank.zeta);
    rec.nOuter = nOuter;
    if strcmp(cfg.loop,'double')
        innerDone = rec.relSd <= cfg.innerTol;
    end
    rec.costGap = (rec.cost - Cstar)/Cstar;
    rec.dMove = NaN; if ~isempty(dPrev), rec.dMove = max(abs(dk - dPrev)); end
    rec.tFit = tFit; rec.tPop = tPop; rec.tOpt = tOpt; rec.tEv = tEv; rec.tAcq = 0;
    rec.elapsed = toc(tStart);
    dPrev = dk;
    if nHF >= cfg.budget
        hist = appendrec(hist,rec); break;
    end
    % ---------- acquisition
    q = min(cfg.B, cfg.budget - nHF);
    t1 = tic;
    aopt = struct('seed',cfg.seed*1000+it);
    switch cfg.method
        case {'stmsur','stmsur_mc','pointsur','stmsur_1d','tmp'}
            if any(strcmp(cfg.method,{'stmsur_1d','tmp'}))
                evs = {ev0}; dw = 1;
            else
                [designs,dw] = neighbourhood(pop,bank,sol,cfg);
                evs = {ev0};
                for jd = 2:size(designs,1)
                    P = c36_pointfield(pop,bank,designs(jd,:));
                    evs{jd} = c36_events(model,bank,P,struct('K',cfg.K,'Emax',cfg.Emax)); %#ok<AGROW>
                end
            end
            mode = struct('stmsur','closed','stmsur_1d','closed','stmsur_mc','mc','pointsur','point','tmp','tmp');
            [Snew,ainfo] = c36_acq_sur(model,bank,evs,dw,q,mode.(cfg.method),aopt);
        case {'stmsur_rb','stmsur_rb1d'}
            % Rao-Blackwellised event SUR (common mode integrated exactly)
            eopt = struct('K',cfg.K,'Emax',cfg.Emax);
            if strcmp(cfg.method,'stmsur_rb1d')
                evs = {c36_events_rb(model,bank,P0,eopt)}; dw = 1;
            else
                [designs,dw] = neighbourhood(pop,bank,sol,cfg);
                evs = cell(1,size(designs,1));
                for jd = 1:size(designs,1)
                    P = c36_pointfield(pop,bank,designs(jd,:));
                    evs{jd} = c36_events_rb(model,bank,P,eopt);
                end
            end
            [Snew,ainfo] = c36_acq_sur(model,bank,evs,dw,q,'closed',aopt);
        case {'rvsur','rvsur_1d'}
            if strcmp(cfg.method,'rvsur_1d')
                funs = {F0}; dw = 1;
            else
                [designs,dw] = neighbourhood(pop,bank,sol,cfg);
                funs = {F0};
                for jd = 2:size(designs,1)
                    P = c36_pointfield(pop,bank,designs(jd,:));
                    funs{jd} = c36_rvfun(model,bank,P,struct('K',cfg.K,'Tmax',cfg.Tmax)); %#ok<AGROW>
                end
            end
            [Snew,ainfo] = c36_acq_rv(model,bank,funs,dw,q,aopt);
        case {'ukb','silk','efflp','pbalc','random','epak','slabalq'}
            [Snew,ainfo] = c36_acq_point(model,bank,ev0,q,cfg.method,aopt);
        otherwise
            error('c36_run:method','unknown method %s',cfg.method);
    end
    rec.tAcq = toc(t1);
    rec.acq = ainfo;
    hist = appendrec(hist,rec);
    % ---------- HF evaluations of the batch (one call per snapshot)
    Ynew = c36_loadeffect(Snew,c.pGrid,'hf');
    S = [S; Snew]; Y = [Y, Ynew]; nHF = size(S,1);
    if cfg.verbose
        fprintf('[%s B=%d s=%d] it=%2d nHF=%3d d=(%.4f,%.4f) C=%.4f Pf_true=%.4f%% relErr*=%+.3f eps=%.3f t=%.1fs\n', ...
            cfg.method,cfg.B,cfg.seed,it,nHF-q,dk,sum(dk),100*rec.pfTrue,rec.pfRelErrStar,rec.epsHat,toc(tStart));
    end
end
if model.lfOnly, th = {}; else, th = cellfun(@(x) x.theta,model.comp,'UniformOutput',false); end
res = struct('cfg',cfg,'hist',hist,'S',S,'Cstar',Cstar,'dstar',dstar, ...
    'final',hist(end),'elapsed',toc(tStart),'modelTheta',{th});
if ~isempty(cfg.resDir)
    if ~isfolder(cfg.resDir), mkdir(cfg.resDir); end
    f = fullfile(cfg.resDir,sprintf('%s.mat',cfg.label));
    save(f,'res');
end
end

function [designs,dw] = neighbourhood(pop,bank,sol,cfg)
% design neighbourhood on the frozen surrogate constraint boundary:
% incumbent (weight 1/2) and the boundary points at h -/+ dh (1/4 each)
c = c36_config(); a = c.a;
designs = sol.d; dw = 0.5;
for dh = [-cfg.dh cfg.dh]
    h = sol.h + dh;
    q = c36_margin(pop.E,c.pGrid,[h -h],bank.zeta);
    f = @(m) mean(0.5*erfc((m+q)/(sqrt(2)*a))) - c.target;
    mlo = max(0.4-h,0.4+h); mhi = min(1.8-h,1.8+h);
    if mlo>mhi || f(mhi)>0, continue; end
    if f(mlo)<=0, m = mlo; else, m = fzero(f,[mlo mhi]); end
    designs(end+1,:) = [m+h, m-h]; %#ok<AGROW>
    dw(end+1) = 0.25; %#ok<AGROW>
end
dw = dw/sum(dw);
end

function hist = appendrec(hist,rec)
if ~isfield(rec,'acq'), rec.acq = struct(); end
if isempty(hist), hist = rec; else, hist(end+1) = orderfields(rec,hist(1)); end
end

function cfg = defaults(cfg)
d = struct('method','stmsur','B',4,'seed',1,'budget',100,'N0',12, ...
    'spec',struct('spatial','op','mean','lf1'),'K',6,'Emax',400,'Tmax',600,'dh',0.04, ...
    'resDir','','truthDir',fullfile(fileparts(mfilename('fullpath')),'..','results','truth'), ...
    'verbose',true,'label','','loop','single','innerTol',0.05,'maxInner',Inf,'evalVal',false);
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(cfg,f{k}), cfg.(f{k}) = d.(f{k}); end
end
if isempty(cfg.label)
    ns = ''; if isfield(cfg.spec,'nonstat') && cfg.spec.nonstat, ns = 'NS'; end
    lp = ''; if strcmp(cfg.loop,'double'), lp = '_DL'; end
    if strcmp(cfg.loop,'double') && isfinite(cfg.maxInner), lp = sprintf('_DLm%d',cfg.maxInner); end
    cfg.label = sprintf('%s_%s%s%s_B%d_s%03d%s',cfg.method,cfg.spec.spatial,cfg.spec.mean,ns,cfg.B,cfg.seed,lp);
end
end
