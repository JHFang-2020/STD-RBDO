function ev = b36_events(gp,E,opts)
%B36_EVENTS  Trajectory failure events at one design (from b36_popeval).
%   Binary events {M_n <= 0}: in Case 2 every aleatory input is part of the
%   history xi_n, so the event is deterministic given the HF field and no
%   aleatory mode is sampled (there is no exactly integrable common mode;
%   Rao-Blackwellisation is therefore not available).
%   Ambiguous trajectories: [mu - kappa*sd, mu + kappa*sd] contains 0; capped
%   at Emax by pi(1-pi).  Returns their window rows, Clark weights, moments,
%   and the posterior s.d. of the working risk
%       Var(Pf) = N^-2 sum_{n,m} [Phi2(alpha_n,alpha_m;rho_nm) - pi_n pi_m],
%   with rho_nm from the (linearised) joint posterior of the Clark margins.
%   ev.cov(i) rows: margin M_n = sum_k w_nk g_nk (row-major: event by event).
if nargin<3, opts = struct(); end
Emax = getf(opts,'Emax',250); kap = getf(opts,'kappa',3); N = getf(opts,'N',4096);
sd = sqrt(max(E.v,0));
amb = find(E.mu - kap*sd < 0 & E.mu + kap*sd > 0 & sd > 1e-9);
pp = E.pi(amb).*(1-E.pi(amb));
[~,o] = sort(pp,'descend'); amb = amb(o(1:min(Emax,numel(amb))));
nE = numel(amb); K = E.K;
rows = reshape(((amb(:)-1)*K + (1:K))',[],1);          % row-major: event by event
ev = struct('d',E.d,'n',E.idx(amb),'loc',amb,'Z',E.Z(rows,:),'w',E.w(amb,:),'mu',E.mu(amb),'v',E.v(amb), ...
    'pi',E.pi(amb),'K',K,'nE',nE);
% ---- posterior s.d. of the working risk at this design
sdPf = 0;
if nE > 0 && ~isempty(gp)
    A = kron(speye(nE),ones(1,K)).*reshape(-E.w(amb,:)',1,[]);   % nE-by-nE*K aggregation (g = c - delta)
    [~,~,Wh] = b36_gp_pred(gp,ev.Z);
    Cd = b36_gp_k(gp,ev.Z,ev.Z) - Wh'*Wh;
    Cm = full(A*Cd*A'); Cm = (Cm+Cm')/2;
    s = sqrt(max(diag(Cm),1e-300));
    R = Cm./(s*s'); R = min(max(R,-1),1);
    al = -ev.mu./sqrt(max(ev.v,1e-300));
    [I,J] = find(triu(true(nE),1));
    cv = b36_phi2cov(al(I),al(J),R(sub2ind([nE nE],I,J)));
    varSum = sum(ev.pi.*(1-ev.pi)) + 2*sum(cv);
    sdPf = sqrt(max(varSum,0))/N;
    ev.Cm = Cm;
end
ev.sdPf = sdPf;
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
