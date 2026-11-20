function res = b36_run(cfg)
%B36_RUN  One single-loop STD-RBDO run of Case 2 (one method, B, seed).
%   Loop:  fit surrogate -> screen the population at the current design ->
%          solve the frozen-surrogate RBDO -> record checkpoint (truth
%          scored outside learning) -> select a batch of HF actions
%          (history, design) -> evaluate (HF) -> append.
%   Identical for every method; only the acquisition (cfg.method), the
%   surrogate spec and (cfg.loop) the coupling differ.
c = b36_config();
cfg = defaults(cfg);
TR = b36_truth_ref(cfg.truthDir);
starSet = TR.T.scr;                                  % LF (Efield) screening at d*_W, theta = 0.25: independent of truth and method
P = b36_pop('work'); P.ext = []; P.dFull = [];
spec = cfg.spec; lfo = struct('variant',getf(spec,'lf','Efield'));
meanKind = getf(spec,'mean','lf');
eopt = struct('lf',lfo,'meanKind',meanKind);
% ---------- initial HF (identical for all methods of a seed)
rs0 = RandStream('mt19937ar','Seed',2000+cfg.seed);
D0 = c.lower + (c.upper-c.lower).*lhs(rs0,cfg.N0,3);
H0 = rs0.randperm(P.N,cfg.N0)';
Z = []; y = []; acts = zeros(0,5);
for i = 1:cfg.N0, [P,Z,y,acts] = hfAction(P,D0(i,:),H0(i),Z,y,acts,lfo,meanKind,0); end
nHF = cfg.N0;
gp = []; hist = struct([]); it = 0; tStart = tic;
dScr = cfg.dStart; dPrev = []; dCur = []; innerDone = false; nInner = 0; nOuter = 0;
lastTruth = []; tScore = 0; S = [];
while true
    it = it + 1;
    t1 = tic; gp = b36_gp_fit(Z,y,spec,gp); tFit = toc(t1);
    if ~isempty(cfg.debugDir), save(fullfile(cfg.debugDir,sprintf('it%02d.mat',it)),'gp','Z','y','acts','dScr','dPrev'); end
    t1 = tic;
    if isempty(S) || any(S.d ~= dScr), [P,S] = screenAt(P,dScr,lfo); end
    if isempty(dPrev), starts = [dScr; cfg.extraStarts];
    else, starts = dPrev; if mod(it,4)==0, starts = [starts; cfg.extraStarts(1,:)]; end
    end
    [sol,P] = b36_solve(gp,P,S,starts,eopt);
    % re-screen at the solution and re-solve until the screening design and
    % the solution agree (the frozen moving set is only valid near S.d)
    for pass = 1:3
        dS = S.d; [P,S] = screenAt(P,sol.d,lfo);
        if max(abs(sol.d - dS)) < 0.003 && all(ismember(S.act,sol.act0)), break; end
        st = sol.d; if ~isempty(dPrev), st = [st; dPrev]; end %#ok<AGROW>
        [sol,P] = b36_solve(gp,P,S,st,eopt);
    end
    if max(abs(sol.d - S.d)) > 0, [P,S] = screenAt(P,sol.d,lfo); end
    tOpt = toc(t1);
    if cfg.verbose && ~isempty(sol.msgs), fprintf('  solver messages: %s\n',strjoin(unique(sol.msgs),' | ')); end
    if cfg.verbose && ~sol.feasible, fprintf('  WARNING: no feasible surrogate design found (pf %.4f)\n',sol.pf); end
    dk = sol.d;
    eo = eopt; eo.keepAll = true;
    if strcmp(cfg.loop,'double')
        % nested coupling with data reuse: the reliability analysis stays at
        % the current design dCur until its risk s.d. <= innerTol*target (or
        % maxInner batches); then the design moves to the surrogate optimum
        if isempty(dCur), dCur = sol.d; nOuter = 1; nInner = 0; end
        dk = dCur;
        if max(abs(dk - S.d)) > 0, [P,S] = screenAt(P,dk,lfo); end
    end
    % ---------- events at the incumbent (shared diagnostics for all methods)
    t1 = tic;
    [E0,P] = b36_popeval(gp,P,dk,S.act,eo);
    ev0 = b36_events(gp,E0,struct('N',P.N));
    if strcmp(cfg.loop,'double') && (ev0.sdPf/c.target <= cfg.innerTol || nInner >= cfg.maxInner) && max(abs(sol.d - dCur)) > 0
        dCur = sol.d; dk = dCur; nOuter = nOuter + 1; nInner = 0;      % inner analysis converged: move now
        [P,S] = screenAt(P,dk,lfo);
        [E0,P] = b36_popeval(gp,P,dk,S.act,eo);
        ev0 = b36_events(gp,E0,struct('N',P.N));
    end
    tEv = toc(t1);
    % ---------- checkpoint metrics (truth never fed back)
    t1 = tic;
    if isempty(lastTruth) || max(abs(lastTruth.d - dk)) > 0
        % the truth screening always uses the same (Efield) simplified physics,
        % independent of the method's surrogate specification
        if strcmp(lfo.variant,'Efield'), St = S; else, [P,St] = b36_screen(P,dk,c.thetaActive,struct()); end
        lastTruth = b36_truth_at(P,St,struct('theta',0.12));
        if ~lastTruth.exact, lastTruth = b36_truth_at(P,St,struct('theta',0.25,'maxHF',2000)); end
    end
    T = lastTruth;
    [Es,P] = b36_popeval(gp,P,TR.dW,starSet,eopt);
    [~,ic] = ismember(TR.crit,Es.idx);
    tTruth = toc(t1);
    rec = struct();
    rec.iter = it; rec.nHF = nHF; rec.d = dk; rec.cost = sum(dk);
    rec.costGapW = (rec.cost - TR.CW)/TR.CW;
    rec.qTrue = T.q; rec.qErr = T.qErr; rec.nFailW = T.nFail; rec.pfTrueW = T.pf;
    rec.feasErrW = (T.pf - c.target)/c.target; rec.truthExact = T.exact; rec.nTruthHF = T.nHF;
    rec.pfSur = E0.pf; rec.pfPlug = sum(E0.hm<=0)/P.N; rec.sdPf = ev0.sdPf; rec.relSd = ev0.sdPf/c.target;
    rec.pfSurStar = sum(Es.pi)/P.N;                  % starSet: LF screening at d* (method- and truth-independent)
    rec.pfRelErrStar = (rec.pfSurStar - TR.pfW)/TR.pfW;
    rec.hRMSEcrit = rms(Es.hm(ic) - TR.margCrit);    % plug-in (posterior-mean) margins of the near-limit histories
    rec.misclassStar = sum((Es.hm(ic)<=0) ~= (TR.margCrit<=0));
    rec.nAmbiguous = ev0.nE; rec.nActive = numel(S.act); rec.solFeasible = sol.feasible; rec.solErrors = numel(sol.msgs);
    rec.nOuter = nOuter; rec.dMove = NaN; if ~isempty(dPrev), rec.dMove = max(abs(dk - dPrev)); end
    rec.gpTheta = gp.theta; rec.gpNu = gp.nu; rec.gpGamma = gp.gamma;
    rec.tFit = tFit; rec.tOpt = tOpt; rec.tEv = tEv; rec.tTruth = tTruth; rec.tAcq = 0;
    tScore = tScore + tTruth;
    rec.elapsed = toc(tStart) - tScore;               % method wall time (truth scoring excluded)
    dPrev = dk; dScr = dk;
    if nHF >= cfg.budget
        hist = appendrec(hist,rec); break;
    end
    % ---------- acquisition (exact repeats of an HF action are never candidates)
    q = min(cfg.B, cfg.budget - nHF);
    t1 = tic; aseed = cfg.seed*1000 + it;
    rx = RandStream('mt19937ar','Seed',aseed);
    doneAt = @(D,H) ismember([D H],acts(:,1:4),'rows');
    switch cfg.method
        case {'stmsur','stmsur_1d'}
            if strcmp(cfg.method,'stmsur_1d') || strcmp(cfg.loop,'double')
                designs = dk; dw = 1;
            else
                [designs,dw,P] = neighbourhood(gp,P,S,dk,cfg.dh,eopt);
            end
            evs = cell(1,size(designs,1)); Cs = cell(1,size(designs,1)+1);
            for j = 1:size(designs,1)
                if j==1, Ej = E0; evj = ev0;
                else, [Ej,P] = b36_popeval(gp,P,designs(j,:),S.act,eo); evj = b36_events(gp,Ej,struct('N',P.N)); end
                evs{j} = evj;
                if evj.nE > 0
                    [~,o] = sort(evj.pi.*(1-evj.pi),'descend');
                    Cs{j} = b36_cands(Ej,P,evj.loc(o(1:min(cfg.Ccap,evj.nE))));
                end
            end
            Cs{end} = b36_cands(E0,P,rx.randperm(numel(E0.idx),min(cfg.nExplore,numel(E0.idx)))');
            Cand = b36_catcands(Cs);
            Cand = b36_subcands(Cand,~doneAt(Cand.d,Cand.hist));
            [selC,ainfo] = b36_acq_sur(gp,evs,dw,Cand,q);
            newD = Cand.d(selC,:); newH = Cand.hist(selC);
        case {'ukb','silk','slabalq','epak','pbalc','efflp','random'}
            excl = doneAt(repmat(dk,numel(E0.idx),1),E0.idx);
            [sel,ainfo] = b36_acq_point(gp,E0,P,q,cfg.method,aseed,excl);
            newD = repmat(dk,numel(sel),1); newH = E0.idx(sel);
        otherwise
            error('b36_run:method','unknown method %s',cfg.method);
    end
    if numel(newH) < q                % nothing informative left: fill with random new actions at dk
        free = E0.idx(~doneAt(repmat(dk,numel(E0.idx),1),E0.idx) & ~ismember(E0.idx,newH(ismember(newD,dk,'rows'))));
        add = free(rx.randperm(numel(free),min(q-numel(newH),numel(free))));
        newD = [newD; repmat(dk,numel(add),1)]; newH = [newH; add]; ainfo.filled = numel(add);
    end
    rec.tAcq = toc(t1); rec.acq = ainfo;
    hist = appendrec(hist,rec);
    % ---------- HF evaluations of the batch (one call per action)
    for i = 1:numel(newH)
        [P,Z,y,acts] = hfAction(P,newD(i,:),newH(i),Z,y,acts,lfo,meanKind,it);
    end
    if strcmp(cfg.loop,'double'), nInner = nInner + 1; end
    nHF = nHF + numel(newH);
    if cfg.verbose
        fprintf('[%s B=%d s=%d] it=%2d nHF=%3d d=(%.4f %.4f %.4f) C=%.5f gap=%+.4f%% qErr=%+.3f%% nFail=%d pfSur=%.4f sd=%.4f t=%.0fs\n', ...
            cfg.method,cfg.B,cfg.seed,it,rec.nHF,dk,rec.cost,100*rec.costGapW,100*rec.qErr,rec.nFailW,rec.pfSur,rec.relSd,toc(tStart));
    end
end
res = struct('cfg',cfg,'hist',hist,'acts',acts,'final',hist(end),'elapsed',toc(tStart),'gpTheta',gp.theta, ...
    'TRd',TR.dW,'TRC',TR.CW);
if ~isempty(cfg.resDir)
    if ~isfolder(cfg.resDir), mkdir(cfg.resDir); end
    save(fullfile(cfg.resDir,sprintf('%s.mat',cfg.label)),'res');
end
end

function [P,Z,y,acts] = hfAction(P,d,h,Z,y,acts,lfo,meanKind,it)
% one HF call: history h at design d -> its block of Kobs windows
[P.cache,ld] = b36_cache(P.cache,P.N,P.xi,h);
blk = b36_block(d,P.xi(h,:),ld,lfo);
o = b36_hf(d,P.xi(h,:));
if strcmp(meanKind,'qs'), Wm = blk.Wqs; Zb = blk.Z; Zb(:,end) = 0; else, Wm = blk.Wlf; Zb = blk.Z; end
Z = [Z; Zb]; y = [y; log(o.W(blk.wk)'./Wm')];
acts = [acts; d h it];
end

function [designs,dw,P] = neighbourhood(gp,P,S,dk,dh,eopt)
% design neighbourhood on the frozen surrogate constraint boundary: the
% incumbent (weight 1/2) and the boundary points obtained by moving d1 or d2
% by +/- dh (inside the box) and restoring Pf = target with d3 (monotone root)
c = b36_config();
designs = dk; dw = 0.5;
eopt.sFloor = 1e-3;                        % same smoothing as the design solver
pfd = @(x) b36_popeval(gp,P,x,S.act,eopt);
for i = 1:2
    for s = [-1 1]
        d = dk; d(i) = d(i) + s*dh;
        if d(i) < c.lower(i)-1e-12 || d(i) > c.upper(i)+1e-12, continue; end
        f = @(x3) getpf(pfd([d(1:2) x3])) - c.target;
        lo = max(c.lower(3),dk(3)-0.03); hi = min(c.upper(3),dk(3)+0.03);
        flo = f(lo); fhi = f(hi);
        if fhi > 0 || flo < 0, continue; end
        d(3) = fzero(f,[lo hi],optimset('TolX',1e-5));
        designs(end+1,:) = d; %#ok<AGROW>
    end
end
if size(designs,1) > 1, dw = [0.5, 0.5*ones(1,size(designs,1)-1)/(size(designs,1)-1)]; else, dw = 1; end
end

function p = getpf(E), p = E.pf; end

function [P,S] = screenAt(P,d,lfo)
% full screening after a design move > c.dFullMove, else the extended set only
c = b36_config();
if isempty(P.dFull) || max(abs(d - P.dFull)) > c.dFullMove
    [P,S] = b36_screen(P,d,c.thetaActive,lfo); P.ext = S.ext; P.dFull = d;
else
    [P,S] = b36_screen(P,d,c.thetaActive,lfo,P.ext);
end
end

function U = lhs(rs,n,k)
U = (cell2mat(arrayfun(@(j) rs.randperm(n)',1:k,'UniformOutput',false)) - rs.rand(n,k))/n;
end

function hist = appendrec(hist,rec)
if ~isfield(rec,'acq'), rec.acq = struct(); end
if isempty(hist), hist = rec; else, hist(end+1) = orderfields(rec,hist(1)); end
end

function cfg = defaults(cfg)
d = struct('method','stmsur','B',4,'seed',1,'budget',100,'N0',16, ...
    'spec',struct('mean','lf','nonstat',true,'lf','Efield'),'dh',0.008,'Ccap',120,'nExplore',30, ...
    'dStart',[0.17 0.29 0.20],'extraStarts',[0.16 0.30 0.19; 0.18 0.28 0.20], ...
    'resDir','','truthDir',fullfile(fileparts(mfilename('fullpath')),'..','results','truth'), ...
    'verbose',true,'label','','loop','single','innerTol',0.05,'maxInner',Inf,'debugDir','');
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(cfg,f{k}), cfg.(f{k}) = d.(f{k}); end
end
if isempty(cfg.label)
    sp = [getf(cfg.spec,'lf','Efield') ternary(getf(cfg.spec,'nonstat',false),'NS','')];
    lp = ''; if strcmp(cfg.loop,'double'), lp = '_DL'; if isfinite(cfg.maxInner), lp = sprintf('_DLm%d',cfg.maxInner); end, end
    cfg.label = sprintf('%s_%s_B%d_s%03d%s',cfg.method,sp,cfg.B,cfg.seed,lp);
end
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
