function V = b36_pf_val(d,theta)
%B36_PF_VAL  Failure probability of one design on the validation bank for
%   designs that may be far from the reliability target (deterministic
%   optimum, LF-only designs).  The simplified physics screens all c.Nval
%   histories; HF (full grid) is run for every history whose LF working peak
%   lies within log-margin theta of 160 MPa.  Histories below the band are
%   certified safe and histories above it certified failing when the band's
%   largest |log HF/LF| (+0.03) does not reach across the band edge.
if nargin<2, theta = 0.15; end
c = b36_config(); Pv = b36_pop('val'); N = Pv.N;
pkLF = zeros(N,1);
for k0 = 1:256:N
    ii = (k0:min(N,k0+255))';
    o = b36_lf(d,Pv.xi(ii,:),b36_loads(Pv.xi(ii,:)),struct());
    pkLF(ii) = max(o.W,[],2);
end
lo = c.allow*exp(-theta); hi = c.allow*exp(theta);
band = find(pkLF >= lo & pkLF <= hi); xi = Pv.xi(band,:); pk = zeros(numel(band),1);
parfor i = 1:numel(band)
    h = b36_hf(d,xi(i,:),'full'); pk(i) = h.peak;
end
e = max(abs(log(pk./pkLF(band))));
nAbove = nnz(pkLF > hi);
nFail = nnz(pk >= c.allow) + nAbove; pf = nFail/N;
cert = lo*exp(e + 0.03) < c.allow && (nAbove==0 || hi*exp(-e - 0.03) >= c.allow);
V = table(d(1),d(2),d(3),sum(d),pf,sqrt(pf*(1-pf)/N),nFail,cert,numel(band),nAbove, ...
    'VariableNames',{'d1','d2','d3','cost','pfVal','pfValSE','nFailVal','certified','nHF','nAboveBand'});
end
