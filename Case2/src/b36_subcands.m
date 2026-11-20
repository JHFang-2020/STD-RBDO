function C = b36_subcands(C,keep)
%B36_SUBCANDS  Keep the candidate HF actions with keep(i) = true.
keep = logical(keep(:)); u = find(keep); Ko = C.Ko;
rows = reshape(((u(:)-1)*Ko + (1:Ko))',[],1);
C.d = C.d(u,:); C.hist = C.hist(u); C.wk = C.wk(u,:); C.Wlf = C.Wlf(u,:); C.Wqs = C.Wqs(u,:);
C.Z = C.Z(rows,:); C.nC = numel(u);
end
