function C = b36_catcands(Cs)
%B36_CATCANDS  Concatenate candidate sets (cell array of b36_cands outputs).
Cs = Cs(~cellfun(@isempty,Cs));
C = Cs{1};
for k = 2:numel(Cs)
    C.d = [C.d; Cs{k}.d]; C.hist = [C.hist; Cs{k}.hist]; C.wk = [C.wk; Cs{k}.wk];
    C.Wlf = [C.Wlf; Cs{k}.Wlf]; C.Wqs = [C.Wqs; Cs{k}.Wqs]; C.Z = [C.Z; Cs{k}.Z];
    C.nC = C.nC + Cs{k}.nC;
end
% the same HF action (same history at the same design) only once
[~,u] = unique([C.d C.hist],'rows','stable');
if numel(u) < C.nC
    Ko = C.Ko; rows = reshape(((u(:)-1)*Ko + (1:Ko))',[],1);
    C.d = C.d(u,:); C.hist = C.hist(u); C.wk = C.wk(u,:); C.Wlf = C.Wlf(u,:); C.Wqs = C.Wqs(u,:);
    C.Z = C.Z(rows,:); C.nC = numel(u);
end
end
