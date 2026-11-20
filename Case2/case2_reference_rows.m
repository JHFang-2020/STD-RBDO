function R = case2_reference_rows(nWorkers)
%CASE2_REFERENCE_ROWS  First three rows of Table 5 (no learning involved):
%   - deterministic optimization (b36_ddo): nominal inputs (xi = 0: nominal
%     geometry, modulus and load histories), min d1 + d2 + d3 subject to a
%     nominal HF peak stress (full grid) of at most 160 MPa, SQP within the
%     design bounds [0.16, 0.30]^3 m;
%   - simplified-model-only RBDO (b36_lf_only, variant 'Efield'): d1 and d2
%     at their bounds (0.16 m, 0.30 m), as at the reference solution, and the
%     smallest d3 for which the 41st-largest simplified-model peak stress of
%     the 4096 working histories equals 160 MPa (Pf = 1 % with the simplified
%     model alone; no HF simulation);
%   - reference solution (b36_truth_ref): exact (HF) optimum of the 4096
%     working histories (Appendix A of the paper).
%   All probabilities of failure are evaluated with the HF model (full grid)
%   on the independent validation population (16,384 histories):
%   b36_pf_val for the deterministic and simplified-model designs (far from
%   the target), b36_validate for the reference solution.
%   Everything is cached and computed only when missing (nWorkers parallel
%   workers): results/truth/truth_ref.mat, results/C2_ddo.mat,
%   results/C2_lf_only_designs.csv, results/C2_val_cache.mat.
%   Writes results/reference_rows.mat.
if nargin<1 || isempty(nWorkers), nWorkers = 4; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'src'));
resDir = fullfile(here,'results');
if ~isfolder(resDir), mkdir(resDir); end
% ---- reference solution: exact optimum of the working population
if ~isfile(fullfile(resDir,'truth','truth_ref.mat')), b36_pool(nWorkers); end
TR = b36_truth_ref();
% ---- deterministic optimum
fd = fullfile(resDir,'C2_ddo.mat');
if isfile(fd), q = load(fd); q = q.out; else, b36_pool(nWorkers); q = b36_ddo(); end
% ---- RBDO with the simplified physical model only
fl = fullfile(resDir,'C2_lf_only_designs.csv');
if ~isfile(fl), b36_pool(nWorkers); b36_lf_only({'Efield'}); end
LF = readtable(fl); lf = LF(strcmp(LF.lf,'Efield'),:);
dLF = [lf.d_1 lf.d_2 lf.d_3];
% ---- validation probabilities of failure (cached)
cf = fullfile(resDir,'C2_val_cache.mat');
if isfile(cf), s = load(cf); Vc = s.Vc; else, Vc = emptyCache(); end
if isnan(lookupPf(Vc,TR.dW))
    b36_pool(nWorkers);
    V = b36_validate(TR.dW); Vc = [Vc; V(:,Vc.Properties.VariableNames)]; save(cf,'Vc');
end
if isnan(lookupPf(Vc,dLF))
    b36_pool(nWorkers);
    VL = b36_pf_val(dLF); VL.qVal = NaN; Vc = [Vc; VL(:,Vc.Properties.VariableNames)]; save(cf,'Vc');
end
R = table(["Deterministic optimization";"Simplified-model-only RBDO";"Reference solution (exact model)"], ...
    ["—";"0";"—"],[q.d(1);dLF(1);TR.dW(1)],[q.d(2);dLF(2);TR.dW(2)],[q.d(3);dLF(3);TR.dW(3)], ...
    [q.cost;sum(dLF);TR.CW],[q.pf;lookupPf(Vc,dLF);lookupPf(Vc,TR.dW)], ...
    'VariableNames',{'name','hf','d1','d2','d3','C','pf'});
save(fullfile(resDir,'reference_rows.mat'),'R');
disp(R);
end

function p = lookupPf(Vc,dd)
[tf,loc] = ismembertol(dd,[Vc.d1 Vc.d2 Vc.d3],1e-12,'ByRows',true,'DataScale',1);
p = nan(size(dd,1),1); p(tf) = Vc.pfVal(loc(tf));
end

function Vc = emptyCache()
z = zeros(0,1);
Vc = table(z,z,z,z,z,z,z,z,false(0,1),z,'VariableNames', ...
    {'d1','d2','d3','cost','pfVal','pfValSE','nFailVal','qVal','certified','nHF'});
end
