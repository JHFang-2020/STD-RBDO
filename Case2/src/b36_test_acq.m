function b36_test_acq()
%B36_TEST_ACQ  Unit tests of the Case 2 event / acquisition assembly against
%   brute-force GP conditioning (regression test for ordering defects of the
%   Case 1 type): (1) linearised event-margin covariance, (2) explained
%   variances after a greedy batch of block observations (exact Schur
%   updates), (3) first-pick gain, (4) Kriging-Believer window variances.
c = b36_config();
P = b36_pop('work'); d0 = [0.17 0.29 0.20];
[P,S] = b36_screen(P,d0,c.thetaActive);
rs = RandStream('mt19937ar','Seed',3);
hist = S.act(rs.randperm(numel(S.act),40));
D = d0 + 0.01*(2*rs.rand(40,3)-1);
Z = []; y = [];
for i = 1:40
    [P.cache,ld] = b36_cache(P.cache,P.N,P.xi,hist(i));
    blk = b36_block(D(i,:),P.xi(hist(i),:),ld); o = b36_hf(D(i,:),P.xi(hist(i),:));
    Z = [Z; blk.Z]; y = [y; log(o.W(blk.wk)'./blk.Wlf')]; %#ok<AGROW>
end
gp = b36_gp_fit(Z,y,struct('mean','lf','nonstat',true));
[E,P] = b36_popeval(gp,P,d0,S.act,struct('keepAll',true));
ev = b36_events(gp,E,struct('N',P.N));
fprintf('events %d, sdPf %.3e\n',ev.nE,ev.sdPf);
X = gp.Z; Kxx = bf_k(gp,X,X) + gp.tau2*eye(size(X,1));
% (1) event-margin covariance (linearised, Clark weights), events in order
A = kron(eye(ev.nE),ones(1,ev.K)).*reshape(-ev.w',1,[]);
for e = 1:ev.nE        % rows of ev.Z must be the windows of event e, in order
    r = (e-1)*ev.K + (1:ev.K);
    assert(max(max(abs(ev.Z(r,4:20) - P.xi(ev.n(e),:)))) < 1e-15,'event rows mixed between histories');
    assert(max(abs(ev.Z(r,21)' - (c.tWin(E.wk(ev.loc(e),:))-1)/3)) < 1e-15,'event windows out of order');
end
sub = 1:min(60,ev.nE); rr = reshape(((sub(:)-1)*ev.K + (1:ev.K))',[],1);
Cd = postcov(gp,ev.Z(rr,:),ev.Z(rr,:),X,Kxx);
Cm = A(sub,rr)*Cd*A(sub,rr)';
e1 = max(abs(Cm(:) - reshape(ev.Cm(sub,sub),[],1)))/max(abs(Cm(:)));
fprintf('(1) event covariance rel error %.2e\n',e1); assert(e1 < 1e-8);
% (2) batch SUR: Schur-updated explained variances vs brute force
[~,o] = sort(ev.pi.*(1-ev.pi),'descend');
C = b36_catcands({b36_cands(E,P,ev.loc(o(1:min(30,end)))), b36_cands(E,P,rs.randperm(numel(E.idx),10)')});
[sel,info] = b36_acq_sur(gp,{ev},1,C,3);
Ko = C.Ko; Ysel = []; for b = 1:numel(sel), Ysel = [Ysel; C.Z((sel(b)-1)*Ko + (1:Ko),:)]; end %#ok<AGROW>
Xa = [X; Ysel];
Ka = [Kxx, bf_k(gp,X,Ysel); bf_k(gp,Ysel,X), bf_k(gp,Ysel,Ysel) + gp.tau2*eye(size(Ysel,1))];
allr = (1:ev.nE*ev.K)';
Cb = postcov(gp,ev.Z,ev.Z,X,Kxx); Ca = postcov(gp,ev.Z,ev.Z,Xa,Ka);
vb = diag(A*Cb*A'); va = diag(A*Ca*A');
e2 = max(abs((vb - va) - info.s2E))/max(vb);
fprintf('(2) explained variances after %d picks: max error %.2e (relative to max var)\n',numel(sel),e2); assert(e2 < 1e-8);
% (3) gain of the first pick by brute force
Y1 = C.Z((sel(1)-1)*Ko + (1:Ko),:);
K1 = [Kxx, bf_k(gp,X,Y1); bf_k(gp,Y1,X), bf_k(gp,Y1,Y1) + gp.tau2*eye(Ko)];
v1 = diag(A*postcov(gp,ev.Z,ev.Z,[X; Y1],K1)*A');
al = -ev.mu./sqrt(ev.v);
g1 = sum(b36_phi2gain(al,min((vb-v1)./ev.v,1-1e-12)));
fprintf('(3) first-pick gain %.6e (algorithm) vs %.6e (brute force)\n',info.scores(1),g1);
assert(abs(g1 - info.scores(1)) < 1e-8*max(abs(g1),1e-12));
% (4) Kriging Believer variances of the tracked windows
[selK,ik] = b36_acq_point(gp,E,P,3,'ukb',7);
Ck = b36_cands(E,P,selK); Yk = Ck.Z;
Kk = [Kxx, bf_k(gp,X,Yk); bf_k(gp,Yk,X), bf_k(gp,Yk,Yk) + gp.tau2*eye(size(Yk,1))];
chk = rs.randperm(numel(E.idx),40);
rows = reshape(((chk(:)-1)*E.K + (1:E.K))',[],1);
vbf = diag(postcov(gp,E.Z(rows,:),E.Z(rows,:),[X; Yk],Kk));
vk = reshape(ik.vtFinal(chk,:)',[],1);
e4 = max(abs(vbf - vk))/max(vbf);
fprintf('(4) Kriging-Believer window variances: max rel error %.2e\n',e4); assert(e4 < 1e-8);
disp('b36_test_acq: all checks passed');
end

function C = postcov(gp,A,B,X,Kxx)
C = bf_k(gp,A,B) - bf_k(gp,A,X)/Kxx*bf_k(gp,X,B);
end

function K = bf_k(gp,A,B)
% independent kernel (vectorised over pairs, no shared code with b36_gp_k)
grp = [1 2 3, 4 4 4, 5 5, 6*ones(1,12), 7, 8];
R2 = zeros(size(A,1),size(B,1));
for j = 1:size(A,2)
    R2 = R2 + ((A(:,j) - B(:,j)')/gp.ell(grp(j))).^2;
end
r = sqrt(5*R2);
K = gp.sf2*(1 + r + r.^2/3).*exp(-r).*(sqrt(1 + gp.gamma*A(:,end).^2)*sqrt(1 + gp.gamma*B(:,end).^2)');
end
