function C = b36_cands(E,P,loc)
%B36_CANDS  Candidate HF actions: histories E.idx(loc) at design E.d.
%   Each candidate is one HF call returning its observation block: the
%   Kobs windows with the largest simplified-physics maxima at that design
%   (identical to the block that b36_block builds after selection).
%   C.Z : nC*Kobs-by-22 (row-major: candidate by candidate)
c = b36_config(); Ko = c.Kobs;
loc = loc(:); nC = numel(loc);
[~,ord] = sort(E.WlfAll(loc,:),2,'descend'); wk = ord(:,1:Ko);
lin = sub2ind(size(E.WlfAll),repmat(loc,1,Ko),wk);
Wl = E.WlfAll(lin); Wq = E.WqsAll(lin);
hi = E.idx(loc);
rowsC = repelem((1:nC)',Ko);
Z = b36_feat(E.d,P.xi(hi(rowsC),:),reshape(c.tWin(wk)',[],1),reshape(log(Wl./Wq)',[],1));
C = struct('d',repmat(E.d,nC,1),'hist',hi,'wk',wk,'Wlf',Wl,'Wqs',Wq,'Z',Z,'nC',nC,'Ko',Ko);
end
