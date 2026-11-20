function sol = c36_solve_rbdo(E,pGrid,zeta,target,box,nh)
%C36_SOLVE_RBDO  min d1+d2  s.t.  Pf(d) <= target  on a frozen envelope bank.
%   Parametrisation d = (m0+h, m0-h): cost = 2*m0 and, for fixed h, Pf is
%   strictly decreasing in m0 (m0 shifts every margin uniformly).  Hence
%   the problem is solved globally by a 1-D search over h with an exact
%   monotone root in m0 for every h (no local-optimum or SQP-tolerance issue).
%   box = [lo1 lo2; hi1 hi2] optional design box (e.g. a trust region).
c = c36_config(); a = c.a;
if nargin<4 || isempty(target), target = c.target; end
if nargin<5 || isempty(box), box = [c.designLower; c.designUpper]; end
if nargin<6 || isempty(nh), nh = 57; end
lo = box(1,:); hi = box(2,:);
hmin = (lo(1)-hi(2))/2; hmax = (hi(1)-lo(2))/2;
hmin = max(hmin,-0.7); hmax = min(hmax,0.7);
hg = linspace(hmin,hmax,nh);
m0 = nan(1,nh);
for k = 1:nh
    m0(k) = rootm(hg(k));
end
[best,k] = min(m0);
if ~isfinite(best)
    sol = struct('feasible',false,'d',[NaN NaN],'cost',NaN,'pf',NaN,'h',NaN,'m0',NaN,'hGrid',hg,'m0Grid',m0);
    return;
end
% local refinement of h (the cost curve m0*(h) is smooth)
kl = max(1,k-1); ku = min(nh,k+1);
if ku>kl
    opt = optimset('TolX',1e-6,'Display','off');
    [hb,mb] = fminbnd(@(h) penal(rootm(h)),hg(kl),hg(ku),opt);
    if mb < best, best = mb; hbest = hb; else, hbest = hg(k); end
else
    hbest = hg(k);
end
d = [best+hbest, best-hbest];
pf = c36_pf(E,pGrid,d,zeta);
sol = struct('feasible',true,'d',d,'cost',sum(d),'pf',pf,'h',hbest,'m0',best,'hGrid',hg,'m0Grid',m0);

    function v = penal(m)
        if isfinite(m), v = m; else, v = 1e3; end
    end
    function m = rootm(h)
        q = c36_margin(E,pGrid,[h -h],zeta);           % margins at m0 = 0
        mlo = max(lo(1)-h, lo(2)+h); mhi = min(hi(1)-h, hi(2)+h);
        if mlo > mhi, m = NaN; return; end
        f = @(mm) mean(0.5*erfc((mm+q)/(sqrt(2)*a))) - target;
        if f(mhi) > 0, m = NaN; return; end
        if f(mlo) <= 0, m = mlo; return; end
        m = fzero(f,[mlo mhi],optimset('TolX',1e-10));
    end
end
