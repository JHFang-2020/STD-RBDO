function V = b36_validate(D)
%B36_VALIDATE  Independent validation of designs (rows of D) on the
%   validation bank (c.Nval histories, own seed), with the FULL HF grid
%   (400 x 6001) peak.  Histories are HF-evaluated when their LF working
%   peak is within log-margin 0.15 of 160 MPa (at least k*+40); the others
%   are certified safe when thr*exp(max|log HF/LF| + 0.03) < 160.
%   V : table with pf, its binomial s.e., nFail, q (k*-th largest full peak).
c = b36_config();
nD = size(D,1); pf = nan(nD,1); se = pf; nf = pf; q = pf; ex = false(nD,1); nh = pf;
parfor k = 1:nD
    Pv = b36_pop('val');
    [Pv,Sv] = b36_screen(Pv,D(k,:),0.15);
    T = b36_truth_at(Pv,Sv,struct('theta',0.15,'full',true,'maxHF',4000));
    pf(k) = T.pf; nf(k) = T.nFail; q(k) = T.q; nh(k) = T.nHF;
    se(k) = sqrt(T.pf*(1-T.pf)/Pv.N);
    ex(k) = T.thr*exp(T.maxErr + 0.03) < c.allow;
end
V = table(D(:,1),D(:,2),D(:,3),sum(D,2),pf,se,nf,q,ex,nh,'VariableNames',{'d1','d2','d3','cost','pfVal','pfValSE','nFailVal','qVal','certified','nHF'});
end
