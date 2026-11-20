function blk = b36_block(d,xi,ld,lfopts)
%B36_BLOCK  Observation block of one history at design d (before the HF call).
%   The retained block is the Kobs critical windows with the largest
%   simplified-physics window maxima; their features are known before the
%   HF call, the HF returns the window maxima (b36_hf), and the training
%   target is  y = log(W_HF/W_mean)  at these windows.
%   ld : loads of this history (b36_cache/b36_loads, one row)
if nargin<4, lfopts = struct(); end
c = b36_config();
o = b36_lf(d,xi,ld,lfopts);
[~,ord] = sort(o.W,'descend'); wk = ord(1:c.Kobs);
Z = b36_feat(d,repmat(xi,c.Kobs,1),c.tWin(wk)',log(o.W(wk)./o.Wqs(wk))');
blk = struct('d',d,'xi',xi,'wk',wk,'Z',Z,'Wlf',o.W(wk),'Wqs',o.Wqs(wk),'WlfAll',o.W,'WqsAll',o.Wqs);
end
