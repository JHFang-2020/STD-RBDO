function [sol,P] = b36_solve(gp,P,S,starts,opts)
%B36_SOLVE  min sum(d)  s.t.  Pf(d) <= target  on the FROZEN surrogate.
%   Pf(d) = (1/N) sum_n Phi(-mu_n(d)/sqrt(v_n(d)))  (posterior expected risk,
%   smooth in d).  The histories that can change state are those whose
%   plug-in margin at the screening design S.d lies in [-band, band]; the
%   rest are counted as failed (margin < -band) or safe.  Their critical
%   windows are frozen during the solve so that Pf(d) is continuous.  The
%   caller re-screens at the solution and re-solves if it moved (b36_run).
%   SQP (fmincon) from the given starts.  The stress is monotone decreasing
%   in every d_i (Case 2 landscape check), so the feasible set is upward
%   closed: an SQP end point that is not feasible is repaired by the
%   monotone root along d3 (then along (1,1,1)).  Cheapest feasible end point wins.
if nargin<5, opts = struct(); end
c = b36_config();
band = getf(opts,'band',0.15);
[E0,P] = b36_popeval(gp,P,S.d,S.act,opts);
mov = abs(E0.hm) <= band | (E0.pi > 1e-6 & E0.pi < 1-1e-6);
nConst = sum(E0.hm(~mov) < 0);
sub = S.act(mov); wk = E0.wk(mov,:);
eo = opts; eo.wk = wk; eo.nConst = nConst;
eo.sFloor = getf(opts,'sFloor',1e-3);         % 0.1 % stress: keeps Pf(d) smooth at the finite-difference scale
fo = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1.5e-3, ...
    'MaxIterations',40,'MaxFunctionEvaluations',200,'ConstraintTolerance',1e-4, ...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-8);
best = struct('d',[NaN NaN NaN],'cost',Inf,'pf',NaN,'feasible',false,'flag',NaN,'msg','');
nEval = 0; cacheD = []; cacheF = []; msgs = {};
for s = 1:size(starts,1)
    x0 = min(max(starts(s,:),c.lower),c.upper);
    try
        [x,~,flag] = fmincon(@(x) sum(x),x0,[],[],[],[],c.lower,c.upper,@con,fo);
    catch err
        x = x0; flag = -99; msgs{end+1} = err.message; %#ok<AGROW>
    end
    pfx = pfAt(x);
    if pfx > c.target*(1+1e-3), x = repair(x); pfx = pfAt(x); end
    feas = pfx <= c.target*(1+1e-3);
    cand = struct('d',x,'cost',sum(x),'pf',pfx,'feasible',feas,'flag',flag,'msg','');
    if isnan(best.pf) || (feas && (~best.feasible || cand.cost < best.cost)) || (~best.feasible && ~feas && pfx < best.pf)
        best = cand;
    end
end
sol = best; sol.nEval = nEval; sol.nMoving = numel(sub); sol.nConst = nConst; sol.act0 = S.act; sol.msgs = msgs;

    function x = repair(x)
        % smallest monotone increase restoring Pf <= target (d3 first, then all)
        for u = {[0 0 1], [1 1 1]}
            mv = u{1} > 0 & x < c.upper;              % components that can still increase
            if ~any(mv), continue; end
            uu = u{1}.*mv; smax = min((c.upper(mv) - x(mv))./uu(mv));
            if smax <= 0, continue; end
            u{1} = uu;
            if pfAt(x + smax*u{1}) > c.target, continue; end
            lo = 0; hi = smax;
            for it = 1:30
                mid = 0.5*(lo+hi);
                if pfAt(x + mid*u{1}) <= c.target, hi = mid; else, lo = mid; end
                if hi - lo < 1e-6, break; end
            end
            x = x + hi*u{1}; return;
        end
        x = c.upper;
    end
    function [cc,ceq] = con(x)
        cc = (pfAt(x) - c.target)/c.target; ceq = [];
    end
    function p = pfAt(x)
        if ~isempty(cacheD)
            j = find(all(abs(cacheD - x(:)')<1e-14,2),1);
            if ~isempty(j), p = cacheF(j); return; end
        end
        [Ex,P] = b36_popeval(gp,P,x,sub,eo);
        p = Ex.pf; nEval = nEval + 1;
        cacheD = [cacheD; x(:)']; cacheF = [cacheF; p];
        if size(cacheD,1) > 60, cacheD = cacheD(end-59:end,:); cacheF = cacheF(end-59:end); end
    end
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
