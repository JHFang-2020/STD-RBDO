function bank = c36_bank(N,seed,tGrid)
%C36_BANK  Scrambled-Sobol population of trajectories (zeta, xi, m).
%   column 1  : zeta  (manufacturing difference mode, N(0,1))
%   columns 2-9: xi   (8 load-process coefficients, N(0,1))
%   column 10 : u_m   (common manufacturing mode / a, N(0,1); used only when
%               the common mode is sampled instead of integrated)
%   States are stored time-major for the grid tGrid: S((it-1)*N + n,:).
c = c36_config();
if nargin<3, tGrid = c.tGrid; end
ps = sobolset(10,'Skip',1);
rs = RandStream('mt19937ar','Seed',seed);
prev = RandStream.setGlobalStream(rs);
ps = scramble(ps,'MatousekAffineOwen');   % seed-dependent randomisation
RandStream.setGlobalStream(prev);
U = net(ps,N);
U = min(max(U,1e-12),1-1e-12);
Q = norminv(U);
bank = struct();
bank.N    = N;
bank.seed = seed;
bank.zeta = Q(:,1);
bank.xi   = Q(:,2:9);
bank.um   = Q(:,10);
bank.t    = tGrid(:)';
[Z1,Z2] = c36_process(bank.xi,bank.t);
nt = numel(bank.t);
bank.S = [Z1(:), Z2(:), repelem(bank.t(:),N)];   % (N*nt)-by-3, time-major
bank.nt = nt;
end
