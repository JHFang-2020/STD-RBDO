function ev = c36_events_rb(model,bank,P,opts)
%C36_EVENTS_RB  Rao-Blackwellised trajectory events for ST-MSUR.
%   The design risk is Pf(d) = mean_n q_n with q_n = Phi(-h_n/a), because
%   the common manufacturing mode is integrated exactly.  The event-level
%   quantity whose epistemic uncertainty is reduced is therefore q_n itself
%   (a "soft" trajectory event), not a binary indicator with one sampled m_n.
%   With h_n ~ N(mu_n, v_n) (Clark margin of K critical points):
%       E[q_n] = Phi(alpha_n),        alpha_n = -mu_n/sqrt(a^2+v_n)
%       E_y Var(q_n | y) reduction = Phi2(alpha,alpha;rho) - Phi(alpha)^2,
%       rho = s_n^2(Q)/(a^2+v_n),
%   i.e. the closed form of the binary event SUR with v_n -> a^2 + v_n.
%   The returned struct is drop-in compatible with c36_acq_sur ('closed'):
%   ev.m = mu_n and ev.v = a^2 + v_n.
c = c36_config(); a = c.a;
if nargin<4, opts = struct(); end
K = getf(opts,'K',6); Emax = getf(opts,'Emax',400); kap = getf(opts,'kappa',3);
N = bank.N; nt = bank.nt;
G = P.g;
[gmin,tmin] = min(G,[],2);
lin1 = (1:N)' + N*(tmin-1);
sig1 = c36_pointsigma(model,bank,P,lin1);
tau = max(0.12, 5*prctile(sig1,95)) + 3*a;
cand = find(G < tau | G < gmin + 4*max(sig1,1e-3));
cand = unique([cand; lin1]);
sig = nan(N,nt);
if ~isempty(cand), sig(cand) = c36_pointsigma(model,bank,P,cand); end
sigF = sig; nanMask = isnan(sig);
sig1m = repmat(sig1,1,nt); sigF(nanMask) = sig1m(nanMask);
lcb = G - kap*sigF; ucb = G + kap*sigF;
lower = min(lcb,[],2); upper = min(ucb,[],2);
% epistemically uncertain AND inside the common-mode band
amb = find(lower < kap*a & upper > -kap*a & sig1 > 1e-5);
[~,tl] = min(lcb(amb,:),[],2);
ia = amb + N*(tl-1);
sm = max(sigF(ia),1e-12);
den = sqrt(a^2 + sm.^2);
al = -G(ia)./den;
pm = 0.5*erfc(-al/sqrt(2));
prox = sm./den.*exp(-0.5*al.^2);                % ~ sensitivity of q_n to the epistemic error
[~,o] = sort(prox,'descend');
nE = min(numel(amb),Emax); evN = amb(o(1:nE));
pts = zeros(nE,K); npts = zeros(nE,1);
for e = 1:nE
    n = evN(e); L = lcb(n,:); L(isnan(sig(n,:))) = Inf;
    lm = find(L <= [Inf L(1:end-1)] & L <= [L(2:end) Inf] & isfinite(L));
    [~,oo] = sort(L(lm)); lm = lm(oo);
    sel = lm(1:min(3,numel(lm)));
    if ~isempty(sel), b = sel(1); sel = [sel, b-1, b+1]; end
    [~,oa] = sort(L); sel = [sel, oa(1:min(nt,K+3))];
    sel = sel(sel>=1 & sel<=nt & isfinite(L(min(max(sel,1),nt))));
    sel = unique(sel,'stable'); sel = sel(1:min(K,numel(sel)));
    npts(e) = numel(sel); pts(e,1:npts(e)) = n + N*(sel-1);
end
allLin = reshape(pts.',[],1); allLin = allLin(allLin>0);   % row-wise: event by event
[uLin,~,map] = unique(allLin);
Su = bank.S(uLin,:); ju = P.j(uLin);
[W,~] = c36_cw(model,Su);
mu = cell(nE,1); Cv = cell(nE,1); wgt = cell(nE,1); mE = zeros(nE,1); vE = zeros(nE,1);
pos = 0;
for e = 1:nE
    k = npts(e); ii = map(pos+(1:k)); pos = pos + k;
    C = zeros(k);
    for l = 1:model.r
        cm = model.comp{l};
        cl = c36_ckern(cm,Su(ii,:),Su(ii,:)) - W{l}(:,ii)'*W{l}(:,ii);
        ph = model.Phi(ju(ii),l); C = C + (ph*ph').*cl;
    end
    C = (C+C')/2;
    mu{e} = G(pts(e,1:k))'; Cv{e} = C;
    [mE(e),vE(e),wgt{e}] = c36_clark_min(mu{e},C);
end
piN = 0.5*erfc(gmin./(sqrt(2)*sqrt(a^2 + sig1.^2)));
piN(evN) = 0.5*erfc(mE./sqrt(2*(a^2 + vE)));
ev = struct('Pj',P.j,'pi',piN,'amb',amb,'n',evN,'pts',pts,'npts',npts,'uLin',uLin, ...
    'mu',{mu},'C',{Cv},'w',{wgt},'m',mE,'v',a^2 + vE,'vEpist',vE, ...
    'G',G,'sig',sig,'tau',tau,'nCut',numel(amb)-nE,'rb',true);
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
