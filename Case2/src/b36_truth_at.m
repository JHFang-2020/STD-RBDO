function T = b36_truth_at(P,S,opts)
%B36_TRUTH_AT  Exact (HF) risk of the population at the screened design S.d.
%   Never used inside learning.  HF is run for every history whose LF peak
%   is within log-margin theta of min(160, k*-th largest LF peak), at least
%   k*+40 histories; the others are certified below the order statistic
%   when  thr*exp(max observed |log HF/LF| + 0.02) < q  (flag T.exact).
%   q    : k*-th largest HF working peak, k* = floor(target*N)+1, i.e.
%          the design satisfies Pf <= target  iff  q < 160 MPa
%   nFail, pf : HF failures / N.   peaks : HF working peaks of T.scr.
%   opts.full = true additionally returns the full-grid HF peaks (validation).
if nargin<3, opts = struct(); end
c = b36_config(); N = P.N;
theta = getf(opts,'theta',0.12); maxHF = getf(opts,'maxHF',900);
kst = floor(c.target*N) + 1;
pkLF = max(S.Wlf,[],2);
[srt,o] = sort(pkLF,'descend');
thr = min(c.allow,srt(kst))*exp(-theta);
nScr = max(sum(pkLF >= thr),kst+40); nScr = min(nScr,maxHF);
scr = o(1:nScr); thrUsed = srt(nScr);
d = S.d; xi = P.xi(scr,:);
pk = zeros(nScr,1); pkF = pk; full = getf(opts,'full',false);
if full, hmode = 'full'; else, hmode = 'work'; end
parfor i = 1:nScr
    h = b36_hf(d,xi(i,:),hmode);
    pk(i) = h.peakW; pkF(i) = h.peak;
end
if full, val = pkF; else, val = pk; end
sv = sort(val,'descend');
q = sv(min(kst,nScr));
e = log(pk./pkLF(scr));
T = struct('d',d,'scr',scr,'peaks',pk,'peaksFull',pkF,'q',q,'qErr',q/c.allow-1, ...
    'nFail',sum(val>=c.allow),'pf',sum(val>=c.allow)/N,'kstar',kst,'thr',thrUsed, ...
    'maxErr',max(abs(e)),'exact',thrUsed*exp(max(abs(e))+0.02) < min(q,c.allow),'nHF',nScr);
% exact: no history outside the screened set can exceed q (order statistic)
% nor 160 MPa (failure count)
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
