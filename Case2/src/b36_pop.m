function P = b36_pop(kind)
%B36_POP  Population of load histories (common random numbers).
%   kind 'work' : working population (c.Nwork), shared by all methods/seeds
%        'val'  : independent validation bank (c.Nval)
%   P.xi : N-by-17; P.cache : load cache (b36_cache) for histories that
%   became active (their exact quasi-static moments and modal forces).
c = b36_config();
switch kind
    case 'work', N = c.Nwork; sd = c.workSeed;
    case 'val',  N = c.Nval;  sd = c.valSeed;
end
rs = RandStream('mt19937ar','Seed',sd);
xi = randn(rs,N,17);
P = struct('kind',kind,'N',N,'xi',xi,'cache',b36_cache([],N));
end
