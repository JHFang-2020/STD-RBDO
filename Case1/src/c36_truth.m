function T = c36_truth(kind,seed,resDir)
%C36_TRUTH  Cached benchmark quantities (never used inside learning).
%   kind 'validation' : independent 2^17 Sobol bank, exact envelopes on the
%                       fine grid (401 t x 201 p) + reference optimum d*
%   kind 'work'       : working population of a seed (8192 histories, 41 t)
%                       + exact envelope on the SAME working grid, used to
%                       isolate surrogate error from discretisation error
%   kind 'cert'       : second independent 2^18 bank for final certification
c = c36_config();
if nargin<3, resDir = fullfile(fileparts(mfilename('fullpath')),'..','results','truth'); end
if ~isfolder(resDir), mkdir(resDir); end
switch kind
    case 'validation', f = fullfile(resDir,'truth_validation.mat');
    case 'cert',       f = fullfile(resDir,'truth_certification.mat');
    case 'work',       f = fullfile(resDir,sprintf('truth_work_seed%03d.mat',seed));
end
if isfile(f), T = load(f); T = T.T; return; end
switch kind
    case {'validation','cert'}
        if strcmp(kind,'validation'), N = 2^17; sd = 20261020; else, N = 2^18; sd = 20261021; end
        bank = c36_bank(N,sd,c.tGrid);
        E = single(c36_true_envelope(bank.xi,c.pRef,c.tRef,512));
        T = struct('N',N,'seed',sd,'zeta',bank.zeta,'E',E,'pGrid',c.pRef);
        if strcmp(kind,'validation')
            sol = c36_solve_rbdo(double(E),c.pRef,bank.zeta,c.target,[],81);
            [pf,~,se] = c36_pf(double(E),c.pRef,sol.d,bank.zeta);
            T.dstar = sol.d; T.Cstar = sol.cost; T.pfStar = pf; T.seStar = se; T.sol = sol;
        end
    case 'work'
        bank = c36_bank(8192,seed,c.tGrid);
        np = c.np; N = bank.N; nt = bank.nt;
        F = zeros(np,N,nt,'single');
        for it = 1:nt
            rows = (it-1)*N + (1:N);
            F(:,:,it) = single(c36_loadeffect(bank.S(rows,:),c.pGrid,'hf'));
        end
        T = struct('bank',bank,'Etrue',c36_envelope_field(F));
end
save(f,'T','-v7.3');
end
