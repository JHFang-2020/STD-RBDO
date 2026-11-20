function F = c36_rvfun(model,bank,P,opts)
%C36_RVFUN  Linearised posterior of the working risk Pf(d) at one design.
%   Pf(d) = (1/N) sum_n Phi(-h_n/a),  h_n = min_{t,p} g_n(t,p)  (common mode
%   integrated exactly).  Under the GP posterior, each h_n is represented by
%   a Clark moment-matched minimum of K critical points, h_n ~ N(mu_n, v_n)
%   with tightness weights w_nk, so that (delta method)
%       Pf - E[Pf] ~= sum_n u_n sum_k w_nk (g(x_nk) - E g(x_nk)),
%       u_n = -phi(mu_n/s_n)/(N s_n),  s_n^2 = a^2 + v_n,
%   i.e. Pf is (approximately) a LINEAR functional  L = sum_x c_x D(x)  of the
%   load-effect field with coefficients c_x = u_n w_nk (sign of D absorbed).
%   Its posterior variance, and the variance a snapshot block would explain,
%   then follow exactly from the GP, INCLUDING cross-trajectory covariance.
if nargin<4, opts = struct(); end
K    = getf(opts,'K',6);
Tmax = getf(opts,'Tmax',600);
c = c36_config(); a = c.a; N = bank.N; nt = bank.nt;
G = P.g;                                           % mean margins, no common mode
[hmu,tmin] = min(G,[],2);
lin1 = (1:N)' + N*(tmin-1);
sig1 = c36_pointsigma(model,bank,P,lin1);
s1 = sqrt(a^2 + sig1.^2);
u1 = exp(-0.5*(hmu./s1).^2)./(sqrt(2*pi)*s1)/N;
score = u1.*max(sig1,1e-12);
[~,o] = sort(score,'descend');
keep = o(score(o) > 1e-6*score(o(1)));
T = keep(1:min(Tmax,numel(keep)));
nT = numel(T);
% s.d. at all times of the retained trajectories
linT = T + N*(0:nt-1);                              % nT-by-nt
sigT = reshape(c36_pointsigma(model,bank,P,linT(:)),nT,nt);
lcb = G(T,:) - 3*sigT;
pts = zeros(nT,K); npts = zeros(nT,1);
for e = 1:nT
    L = lcb(e,:);
    lm = find(L <= [Inf L(1:end-1)] & L <= [L(2:end) Inf]);
    [~,oo] = sort(L(lm)); lm = lm(oo);
    sel = lm(1:min(3,numel(lm)));
    if ~isempty(sel), b = sel(1); sel = [sel, b-1, b+1]; end
    [~,oa] = sort(L); sel = [sel, oa(1:min(nt,K+3))];
    sel = sel(sel>=1 & sel<=nt); sel = unique(sel,'stable'); sel = sel(1:min(K,numel(sel)));
    npts(e) = numel(sel); pts(e,1:npts(e)) = T(e) + N*(sel-1);
end
allLin = reshape(pts.',[],1); allLin = allLin(allLin>0);   % row-wise: trajectory by trajectory
[uLin,~,map] = unique(allLin);
Su = bank.S(uLin,:); ju = P.j(uLin);
[W,~] = c36_cw(model,Su);
mE = zeros(nT,1); vE = zeros(nT,1); wk = cell(nT,1); pos = 0; idx = cell(nT,1);
for e = 1:nT
    k = npts(e); ii = map(pos+(1:k)); pos = pos+k; idx{e} = ii;
    C = zeros(k);
    for l = 1:model.r
        cm = model.comp{l};
        cl = c36_ckern(cm,Su(ii,:),Su(ii,:)) - W{l}(:,ii)'*W{l}(:,ii);
        ph = model.Phi(ju(ii),l); C = C + (ph*ph').*cl;
    end
    [mE(e),vE(e),wk{e}] = c36_clark_min(G(pts(e,1:k))',(C+C')/2);
end
sE = sqrt(a^2 + vE);
u = exp(-0.5*(mE./sE).^2)./(sqrt(2*pi)*sE)/N;
% coefficients on unique points (aggregate repeated states/nodes)
coef = zeros(numel(uLin),1);
for e = 1:nT
    coef(idx{e}) = coef(idx{e}) + u(e)*wk{e};
end
% exact posterior variance of the linear functional (sum over components)
varL = 0;
for l = 1:model.r
    cm = model.comp{l};
    al = coef.*model.Phi(ju,l);
    kv = 0;
    for s = 1:2000:numel(uLin)
        ii = s:min(numel(uLin),s+1999);
        kv = kv + al(ii)'*(c36_ckern(cm,Su(ii,:),Su)*al);
    end
    wa = W{l}*al;
    varL = varL + max(kv - wa'*wa,0);
end
F = struct('lin',uLin,'j',ju,'coef',coef,'S',Su,'varL',varL,'sdPf',sqrt(varL), ...
    'sumIndepVar',sum(u.^2.*vE),'T',T,'nT',nT,'meanRiskPost',mean(0.5*erfc(hmu./(sqrt(2)*s1))));
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
