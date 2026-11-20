function R = case1_reference_rows()
%CASE1_REFERENCE_ROWS  First three rows of Table 3 (no learning involved):
%   - deterministic optimization: nominal manufacturing and the nominal load
%     path (Z = 0), linear programme within the design bounds [0.4, 1.8]^2;
%   - simplified-model-only RBDO: the RBDO problem solved with the simplified
%     physical model (Layer 2 alone) on the validation population;
%   - reference solution: exact model on the validation population (Appendix A,
%     one-dimensional search, computed once by c36_truth('validation')).
%   All probabilities of failure are evaluated with the exact model on the
%   independent validation population (2^17 RQMC samples, 401 x 201 grid).
%   Writes results/reference_rows.mat.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
c = c36_config();
V = c36_truth('validation');                         % built and cached on first use
Ev = double(V.E);
% ---- deterministic optimum within the design bounds
p = c.pRef; tt = c.tRef(:);
D = c36_loadeffect([zeros(numel(tt),2) tt],p,'hf');   % nominal load path, np-by-nt
A = -[1-p p]; b = c.offset - max(D,[],2);
opt = optimoptions('linprog','Display','none');
dD = linprog([1 1],A,b,[],[],c.designLower,c.designUpper,opt)';
pfD = c36_pf(Ev,V.pGrid,dD,V.zeta);
% ---- simplified-model-only RBDO
bank = c36_bank(V.N,V.seed,c.tGrid);
N = V.N; np = numel(c.pRef); nt = numel(c.tRef); Elf = zeros(N,np);
for s = 1:256:N
    ii = s:min(N,s+255); n = numel(ii);
    [Z1,Z2] = c36_process(bank.xi(ii,:),c.tRef);
    S = [Z1(:), Z2(:), repelem(c.tRef(:),n)];
    Dl = reshape(c36_loadeffect(S,c.pRef,'lf'),np,n,nt);
    Elf(ii,:) = max(Dl,[],3)';
end
solL = c36_solve_rbdo(Elf,c.pRef,V.zeta,c.target,[],81);
pfL = c36_pf(Ev,V.pGrid,solL.d,V.zeta);
R = table(["Deterministic optimization";"Simplified-model-only RBDO";"Reference solution (exact model)"], ...
    ["—";"0";"—"],[dD(1);solL.d(1);V.dstar(1)],[dD(2);solL.d(2);V.dstar(2)], ...
    [sum(dD);solL.cost;V.Cstar],[pfD;pfL;V.pfStar],'VariableNames',{'name','hf','d1','d2','C','pf'});
save(fullfile(here,'results','reference_rows.mat'),'R');
disp(R);
end
