function S = c36_lhs_states(n,seed)
%C36_LHS_STATES  Maximin Latin-hypercube initial snapshot states in the box.
c = c36_config();
rs = RandStream('mt19937ar','Seed',seed);
best = -Inf; S = [];
for rep = 1:50
    U = (cell2mat(arrayfun(@(k) rs.randperm(n)',1:3,'UniformOutput',false)) - rs.rand(n,3))/n;
    X = c.stateLower + U.*(c.stateUpper - c.stateLower);
    dmin = min(pdist(X./[1 1 4/6]));
    if dmin > best, best = dmin; S = X; end
end
end
