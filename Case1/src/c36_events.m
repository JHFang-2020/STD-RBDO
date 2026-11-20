function ev = c36_events(model,bank,P,opts)
%C36_EVENTS  Trajectory failure events and their Gaussian (Clark) margins.
%   Event of trajectory n:  M_n = m_n + min_t min_p g(n,t,p) <= 0, with the
%   common manufacturing mode m_n = a*u_m,n sampled (aleatory) and the GP
%   posterior providing the epistemic uncertainty.
%   Returns per-trajectory classification probabilities pi (N-by-1) and, for
%   the ambiguous trajectories, K critical points with their joint posterior
%   covariance, moment-matched mean/variance and tightness weights.
c = c36_config();
if nargin<4, opts = struct(); end
K    = getf(opts,'K',6);
Emax = getf(opts,'Emax',400);
kap  = getf(opts,'kappa',3);
N = bank.N; nt = bank.nt;
mN = c.a*bank.um;
G = P.g + mN;                                   % N-by-nt mean margins incl. common mode
% ---- pass 1: s.d. at each trajectory's mean-minimum point
[gmin,tmin] = min(G,[],2);
lin1 = (1:N)' + N*(tmin-1);
sig1 = c36_pointsigma(model,bank,P,lin1);
% ---- pass 2: s.d. at all points that can matter
tau = max(0.12, 5*prctile(sig1,95));
cand = find(G < tau | G < gmin + 4*max(sig1,1e-3) );
cand = unique([cand; lin1]);
sig = nan(N,nt);
if ~isempty(cand)
    sig(cand) = c36_pointsigma(model,bank,P,cand);
end
% pruned points lie far on the safe side; their s.d. is proxied by the
% trajectory's pass-1 value (only used for classification bounds)
sigF = sig; nanMask = isnan(sig);
sig1m = repmat(sig1,1,nt); sigF(nanMask) = sig1m(nanMask);
lcb = G - kap*sigF; ucb = G + kap*sigF;
lower = min(lcb,[],2); upper = min(ucb,[],2);
piN = double(upper < 0);                         % confidently failed
amb = find(lower < 0 & upper >= 0);
% rank ambiguous trajectories by a cheap marginal proxy and cap their number
[~,tl] = min(lcb(amb,:),[],2);
ia = amb + N*(tl-1);
sm = max(sigF(ia),1e-12);
pmarg = 0.5*erfc(G(ia)./(sqrt(2)*sm));
prox = pmarg.*(1-pmarg);
[~,o] = sort(prox,'descend');
nE = min(numel(amb),Emax);
evN = amb(o(1:nE));
% ---- critical points of every retained event
pts = zeros(nE,K); npts = zeros(nE,1);
for e = 1:nE
    n = evN(e);
    L = lcb(n,:);
    L(isnan(sig(n,:))) = Inf;
    lm = find(L <= [Inf L(1:end-1)] & L <= [L(2:end) Inf] & isfinite(L));
    [~,oo] = sort(L(lm)); lm = lm(oo);
    sel = lm(1:min(3,numel(lm)));
    if ~isempty(sel)
        b = sel(1); sel = [sel, b-1, b+1];
    end
    [~,oa] = sort(L); sel = [sel, oa(1:min(nt,K+3))];
    sel = sel(sel>=1 & sel<=nt & isfinite(L(min(max(sel,1),nt))));
    sel = unique(sel,'stable'); sel = sel(1:min(K,numel(sel)));
    npts(e) = numel(sel); pts(e,1:npts(e)) = n + N*(sel-1);
end
% ---- joint posterior covariance of the K points (per event)
allLin = reshape(pts.',[],1); allLin = allLin(allLin>0);   % row-wise: event by event
[uLin,~,map] = unique(allLin);
Su = bank.S(uLin,:); ju = P.j(uLin);
[W,~] = c36_cw(model,Su);
mu = cell(nE,1); Cv = cell(nE,1); wgt = cell(nE,1);
mE = zeros(nE,1); vE = zeros(nE,1);
pos = 0; idxMap = zeros(size(pts));
for e = 1:nE
    k = npts(e); ii = map(pos+(1:k)); pos = pos + k; idxMap(e,1:k) = ii;
    C = zeros(k);
    for l = 1:model.r
        cm = model.comp{l};
        cl = c36_ckern(cm,Su(ii,:),Su(ii,:)) - W{l}(:,ii)'*W{l}(:,ii);
        ph = model.Phi(ju(ii),l);
        C = C + (ph*ph').*cl;
    end
    C = (C+C')/2;
    mu{e} = G(pts(e,1:k))'; Cv{e} = C;
    [mE(e),vE(e),wgt{e}] = c36_clark_min(mu{e},C);
end
piE = 0.5*erfc(mE./sqrt(2*max(vE,1e-300)));
piE(vE<=0) = double(mE(vE<=0)<=0);
piN(evN) = piE;
% trajectories that are ambiguous but not retained: marginal proxy
rest = setdiff(amb,evN);
if ~isempty(rest)
    [~,tl] = min(lcb(rest,:),[],2); ir = rest + N*(tl-1);
    piN(rest) = 0.5*erfc(G(ir)./(sqrt(2)*max(sigF(ir),1e-12)));
end
ev = struct('Pj',P.j,'pi',piN,'amb',amb,'n',evN,'pts',pts,'npts',npts,'uLin',uLin,'idxMap',idxMap, ...
    'mu',{mu},'C',{Cv},'w',{wgt},'m',mE,'v',vE,'piE',piE,'G',G,'sig',sig,'tau',tau,'nCut',numel(amb)-nE);
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
