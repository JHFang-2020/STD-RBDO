function [P,S] = b36_screen(P,d,theta,lfopts,subset)
%B36_SCREEN  Simplified physics of the population at design d.
%   Histories whose LF peak lies within log-margin theta of the allowable
%   stress are "active".  A full screening evaluates all N histories and
%   returns the extended set (log-margin < c.thetaExt), whose loads are
%   kept in the cache; a subset screening (subset given) evaluates only
%   those histories (used while the design stays within c.dFullMove of the
%   last full screening, so no history outside the set can become active).
%   S.Wlf, S.Wqs : N-by-nW window maxima (LF, exact quasi-static) [MPa]
%                  (zero for histories not evaluated)
%   S.act : active histories;  S.ext : extended set;  S.full : full screening
if nargin<4 || isempty(lfopts), lfopts = struct(); end
if nargin<5, subset = []; end
c = b36_config(); N = P.N;
Wlf = zeros(N,c.nW); Wqs = Wlf;
thExt = max(c.thetaExt,theta);
if isempty(subset), rows = (1:N)'; else, rows = subset(:); end
ch = 256;
for k0 = 1:ch:numel(rows)
    ii = rows(k0:min(numel(rows),k0+ch-1));
    have = P.cache.pos(ii)>0;
    if all(have)
        [~,ld] = b36_cache(P.cache,N,P.xi,ii);
    else
        ld = b36_loads(P.xi(ii,:));
    end
    o = b36_lf(d,P.xi(ii,:),ld,lfopts);
    Wlf(ii,:) = o.W; Wqs(ii,:) = o.Wqs;
    a = ii(max(o.W,[],2) >= c.allow*exp(-thExt) & ~have);
    if ~isempty(a)      % store loads of new extended-set histories (no recomputation)
        [~,loc] = ismember(a,ii);
        k1 = P.cache.count + (1:numel(a));
        P.cache.Mqs(:,:,k1) = ld.Mqs(:,:,loc); P.cache.f(k1,:,:) = ld.f(loc,:,:);
        P.cache.pos(a) = k1; P.cache.count = k1(end);
    end
end
pk = max(Wlf,[],2);
act = find(pk >= c.allow*exp(-theta));
if numel(act) < 60                       % very safe design: keep the 60 most loaded histories
    [~,o] = sort(pk,'descend'); act = sort(o(1:60));
    miss = act(P.cache.pos(act)==0);
    if ~isempty(miss), P.cache = b36_cache(P.cache,N,P.xi,miss); end
end
if isempty(subset), ext = find(pk >= c.allow*exp(-thExt)); ext = union(ext,act); else, ext = rows; end
S = struct('d',d,'Wlf',Wlf,'Wqs',Wqs,'act',act,'ext',ext,'theta',theta,'full',isempty(subset));
end
