function [cache,ld] = b36_cache(cache,N,xi,idx)
%B36_CACHE  Load cache of histories (design-independent, b36_loads).
%   cache = b36_cache([],N)            : empty cache for N histories
%   [cache,ld] = b36_cache(cache,N,xi,idx) : loads of rows idx (computed and
%                                        stored when missing), ld rows in idx order
c = b36_config();
if isempty(cache)
    cache = struct('pos',zeros(N,1),'Mqs',zeros(c.nx,c.ntw,0,'single'), ...
        'f',zeros(0,numel(c.tLF),c.J,'single'),'count',0);
    if nargin<3, ld = []; return; end
end
idx = idx(:);
miss = idx(cache.pos(idx)==0); miss = unique(miss,'stable');
if ~isempty(miss)
    L = b36_loads(xi(miss,:));
    k0 = cache.count;
    cache.Mqs(:,:,k0+(1:numel(miss))) = L.Mqs;
    cache.f(k0+(1:numel(miss)),:,:) = L.f;
    cache.pos(miss) = k0 + (1:numel(miss));
    cache.count = k0 + numel(miss);
end
p = cache.pos(idx);
ld = struct('Mqs',cache.Mqs(:,:,p),'f',cache.f(p,:,:),'n',numel(idx));
end
