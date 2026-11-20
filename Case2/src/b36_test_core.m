function b36_test_core()
%B36_TEST_CORE  Unit tests of the Case 2 surrogate core:
%   (1) window maxima of LF/HF on the working grid, (2) GP fit, (3) every
%   within-history posterior covariance and mean of b36_popeval against a
%   brute-force GP posterior written independently, (4) design solve timing.
c = b36_config();
P = b36_pop('work');
d0 = [0.17 0.29 0.20];
t0 = tic; [P,S] = b36_screen(P,d0,c.thetaActive); fprintf('screen N=%d: %.1fs, active %d, cache %d\n',P.N,toc(t0),numel(S.act),P.cache.count);
fprintf('LF Pf at d0: %.4f  (LF peaks >=160: %d)\n',mean(max(S.Wlf,[],2)>=160),sum(max(S.Wlf,[],2)>=160));
% ---- training data: 40 active histories at 3 designs near d0
rs = RandStream('mt19937ar','Seed',3);
hist = S.act(rs.randperm(numel(S.act),min(40,numel(S.act))));
D = d0 + 0.01*(2*rs.rand(numel(hist),3)-1);
Z = []; y = []; t0 = tic; tHF = 0;
for i = 1:numel(hist)
    [P.cache,ld] = b36_cache(P.cache,P.N,P.xi,hist(i));
    blk = b36_block(D(i,:),P.xi(hist(i),:),ld);
    t1 = tic; o = b36_hf(D(i,:),P.xi(hist(i),:)); tHF = tHF + toc(t1);
    Z = [Z; blk.Z]; y = [y; log(o.W(blk.wk)'./blk.Wlf')]; %#ok<AGROW>
    assert(abs(max(o.W) - o.peakW) < 1e-12);
end
fprintf('40 HF blocks: %.1fs (HF %.1fs); y mean %+.4f sd %.4f\n',toc(t0),tHF,mean(y),std(y));
t0 = tic; gp = b36_gp_fit(Z,y,struct('mean','lf','nonstat',true)); fprintf('GP fit n=%d: %.1fs  ell=%s gamma=%.3g nu=%.2e sf=%.4f beta=%+.4f\n',size(Z,1),toc(t0),mat2str(gp.ell,3),gp.gamma,gp.nu,sqrt(gp.sf2),gp.beta);
% ---- popeval and brute force
t0 = tic; [E,P] = b36_popeval(gp,P,d0,S.act); fprintf('popeval %d histories: %.2fs, Pf_mix %.4f, plug-in %.4f\n',numel(S.act),toc(t0),E.pf,mean(E.hm<=0)*numel(S.act)/P.N);
X = gp.Z; n = size(X,1);
Kxx = bf_k(gp,X,X) + gp.tau2*eye(n);
errC = 0; errM = 0;
for i = rs.randperm(numel(S.act),25)
    r = (i-1)*E.K + (1:E.K);
    A = E.Z(r,:);
    Ksx = bf_k(gp,A,X);
    Cbf = bf_k(gp,A,A) - Ksx/Kxx*Ksx';
    mbf = gp.beta + Ksx/Kxx*(gp.y - gp.beta);
    errC = max(errC,max(abs(Cbf(:) - reshape(E.C(:,:,i),[],1)))/max(abs(Cbf(:))));
    m0 = log(c.allow./E.Wlf(i,:))';
    errM = max(errM,max(abs((m0 - mbf) - E.m(i,:)')));
    % the feature rows must belong to history i and its windows (row-major check)
    assert(all(all(abs(A(:,4:20) - P.xi(S.act(i),:)) < 1e-15)));
    assert(all(abs(A(:,21)' - (c.tWin(E.wk(i,:))-1)/3) < 1e-15));
end
fprintf('popeval vs brute force: max rel cov error %.2e, max mean error %.2e\n',errC,errM);
assert(errC < 1e-8 && errM < 1e-8,'popeval posterior mismatch');
t0 = tic; [sol,P] = b36_solve(gp,P,S,[d0; 0.18 0.30 0.18]); fprintf('solve: %.1fs  d=%s cost %.4f pf %.4f feas %d evals %d moving %d const %d\n',toc(t0),mat2str(sol.d,4),sol.cost,sol.pf,sol.feasible,sol.nEval,sol.nMoving,sol.nConst);
disp('b36_test_core: all checks passed');
end

function K = bf_k(gp,A,B)
% independent brute-force kernel: grouped ARD Matern-5/2 with amplitude
grp = [1 2 3, 4 4 4, 5 5, 6*ones(1,12), 7, 8];
K = zeros(size(A,1),size(B,1));
for i = 1:size(A,1)
    for j = 1:size(B,1)
        r2 = 0;
        for g = 1:8, r2 = r2 + sum((A(i,grp==g) - B(j,grp==g)).^2)/gp.ell(g)^2; end
        r = sqrt(5*r2);
        K(i,j) = gp.sf2*(1 + r + r^2/3)*exp(-r)*sqrt(1 + gp.gamma*A(i,end)^2)*sqrt(1 + gp.gamma*B(j,end)^2);
    end
end
end
