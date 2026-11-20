function c36_test_events()
%C36_TEST_EVENTS  Unit test: every event covariance equals the exact posterior
%   covariance of its own critical points (regression test for the
%   point-ordering defect found in review), for binary and RB events and
%   for the linearised risk functional.
c = c36_config();
WRK = c36_truth('work',1); bank = WRK.bank;
S = c36_lhs_states(30,7); Y = c36_loadeffect(S,c.pGrid,'hf');
model = c36_cgp_fit(S,Y,struct('spatial','op','mean','lf1','nonstat',true));
pop = c36_popeval(model,bank);
P = c36_pointfield(pop,bank,[0.55 1.12]);
for kind = {'binary','rb'}
    if strcmp(kind{1},'binary'), ev = c36_events(model,bank,P); else, ev = c36_events_rb(model,bank,P); end
    err = zeros(numel(ev.n),1);
    for e = 1:numel(ev.n)
        k = ev.npts(e); lin = ev.pts(e,1:k)';
        Cx = exactCov(model,bank.S(lin,:),P.j(lin));
        err(e) = max(abs(Cx(:)-ev.C{e}(:)))/max(abs(Cx(:)));
    end
    fprintf('%-6s events: %d, max relative covariance error %.2e\n',kind{1},numel(ev.n),max(err));
    assert(max(err) < 1e-8,'event covariance mismatch');
end
% risk functional: sdPf vs brute force on the same coefficients
F = c36_rvfun(model,bank,P);
Cx = exactCov(model,F.S,F.j);
fprintf('rvfun: sdPf %.4e, brute force %.4e\n',F.sdPf,sqrt(F.coef'*Cx*F.coef));
assert(abs(F.sdPf - sqrt(F.coef'*Cx*F.coef)) < 1e-8*max(F.sdPf,1e-12) + 1e-14);
disp('c36_test_events: all checks passed');
end

function C = exactCov(model,Sx,jx)
n = size(Sx,1); C = zeros(n);
for l = 1:model.r
    cm = model.comp{l};
    W = cm.L\c36_ckern(cm,cm.S,Sx);
    cl = c36_ckern(cm,Sx,Sx) - W'*W;
    ph = model.Phi(jx,l);
    C = C + (ph*ph').*cl;
end
C = (C+C')/2;
end
