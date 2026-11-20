function [E,P] = b36_popeval(gp,P,d,idx,opts)
%B36_POPEVAL  Posterior of the trajectory margins of histories idx at design d.
%   Window margin (log stress):  g_nk = log(160/W_LF,nk) - delta_nk,
%   delta = GP discrepancy (b36_gp_fit) at the K critical windows of the
%   history (largest LF window maxima, or opts.wk fixed by the caller).
%   Trajectory margin M_n = min_k g_nk: Clark moment matching gives
%   mu_n, v_n and tightness weights;  pi_n = P(M_n <= 0) = Phi(-mu_n/sqrt(v_n)).
%   E.pf   : (1/N) sum pi_n over idx  (+ opts.pfConst for fixed histories)
%   E.Z    : n*K-by-22 window features (history-major), E.m : n-by-K means,
%   E.C    : K-by-K-by-n posterior covariances, E.w : n-by-K Clark weights.
%   gp = [] -> simplified physics only (delta = 0, no uncertainty).
if nargin<5, opts = struct(); end
c = b36_config(); K = getf(opts,'K',c.Kev);
n = numel(idx); idx = idx(:);
[P.cache,ld] = b36_cache(P.cache,P.N,P.xi,idx);
o = b36_lf(d,P.xi(idx,:),ld,getf(opts,'lf',struct()));
if isfield(opts,'wk') && ~isempty(opts.wk)
    wk = opts.wk; K = size(wk,2);
else
    [~,o2] = sort(o.W,2,'descend'); wk = o2(:,1:K);
end
lin = sub2ind(size(o.W),repmat((1:n)',1,K),wk);
Wl = o.W(lin); Wq = o.Wqs(lin);
mean0 = getf(opts,'meanKind','lf');
if strcmp(mean0,'qs'), Wm = Wq; else, Wm = Wl; end
rowsH = repelem((1:n)',K);                        % history of every feature row
Z = b36_feat(d,P.xi(idx(rowsH),:),reshape(c.tWin(wk)',[],1),reshape(log(Wl./Wq)',[],1));
if strcmp(mean0,'qs'), Z(:,end) = 0; end
m0 = log(c.allow./Wm);                            % n-by-K
if isempty(gp)
    m = m0; Cv = zeros(K,K,n); mu = min(m,[],2); v = zeros(n,1); w = double(m==mu);
    pin = double(mu<=0);
else
    [md,~,Wh] = b36_gp_pred(gp,Z);
    m = m0 - reshape(md,K,n)';
    % within-history posterior covariances (diagonal blocks only), row-wise:
    % feature row (i-1)*K+k belongs to history i, window k
    Cv = zeros(K,K,n);
    for k1 = 1:K
        r1 = (0:n-1)'*K + k1;
        for k2 = k1:K
            r2 = (0:n-1)'*K + k2;
            ck = b36_gp_kpair(gp,Z(r1,:),Z(r2,:)) - sum(Wh(:,r1).*Wh(:,r2),1)';
            Cv(k1,k2,:) = reshape(ck,1,1,n); Cv(k2,k1,:) = reshape(ck,1,1,n);
        end
    end
    mu = zeros(n,1); v = mu; w = zeros(n,K);
    for i = 1:n
        [mu(i),v(i),w(i,:)] = b36_clark_min(m(i,:)',Cv(:,:,i));
    end
    s0 = getf(opts,'sFloor',0);                   % design solver only: smoothing floor of the margin s.d.
    pin = 0.5*erfc(mu./sqrt(2*max(v + s0^2,1e-300)));
    pin(v + s0^2 <= 0) = double(mu(v + s0^2 <= 0)<=0);
end
pf = (sum(pin) + getf(opts,'nConst',0))/P.N;
E = struct('d',d,'idx',idx,'wk',wk,'Wlf',Wl,'Wqs',Wq,'Z',Z,'m',m,'C',Cv,'mu',mu,'v',v,'w',w, ...
    'pi',pin,'pf',pf,'hm',min(m,[],2),'K',K);
if getf(opts,'keepAll',false), E.WlfAll = o.W; E.WqsAll = o.Wqs; end
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
