function [sel,info] = b36_acq_sur(gp,evs,dw,C,B)
%B36_ACQ_SUR  ST-MSUR for Case 2: closed-form batch event SUR over the
%   design neighbourhood with exact block Schur updates.
%   evs : cell of event sets (b36_events), one per design of the
%         neighbourhood; dw : their weights (sum 1)
%   C   : candidate HF actions (b36_cands/b36_catcands); each candidate is a
%         history at one of the neighbourhood designs and observes its whole
%         block of Kobs correlated window values (with the GP nugget).
%   Event n: M_n = sum_k w_nk g_nk (Clark), g = log(160/W_LF) - delta.
%   Explained variance of a block Y: s_n^2 = c_n' Sigma_Y^{-1} c_n,
%   c_n = Cov(M_n,Y).  Expected reduction of the event-indicator variance:
%       Delta_n = Phi2(alpha_n,alpha_n;rho_n) - Phi(alpha_n)^2,
%       alpha_n = -mu_n/sqrt(v_n),  rho_n = s_n^2/v_n   (Chevalier et al. 2014).
%   Greedy batch; after each pick every covariance is conditioned exactly
%   (rank-Kobs Schur update), so the joint batch value needs no fantasies.
Ko = C.Ko; nC = C.nC; tau2 = gp.tau2;
% ---------------- events (row-major: event by event, K windows each)
Zev = []; ri = []; ci = []; vals = []; mu = []; v = []; EW = []; nE = 0; off = 0;
for j = 1:numel(evs)
    ev = evs{j}; if isempty(ev) || ev.nE==0, continue; end
    K = ev.K;
    Zev = [Zev; ev.Z]; %#ok<AGROW>
    ri = [ri; reshape(repmat(nE + (1:ev.nE),K,1),[],1)]; %#ok<AGROW>   (column also for nE = 1)
    ci = [ci; off + (1:ev.nE*K)']; %#ok<AGROW>
    vals = [vals; reshape(-ev.w',[],1)]; %#ok<AGROW>   % g = c - delta
    mu = [mu; ev.mu]; v = [v; ev.v]; EW = [EW; dw(j)*ones(ev.nE,1)]; %#ok<AGROW>
    nE = nE + ev.nE; off = off + ev.nE*K;
end
info = struct('nEvents',nE,'nCandidates',nC);
sel = zeros(B,1); score = zeros(B,1);
if nE == 0, sel = []; info.fallback = 'no ambiguous event'; return; end
A = sparse(ri,ci,vals,nE,off);
v = max(v,1e-300); al = -mu./sqrt(v);
[~,~,Wev] = b36_gp_pred(gp,Zev);
[~,~,Wc]  = b36_gp_pred(gp,C.Z);
CE = zeros(nE,nC*Ko);
ch = 1200;
for c0 = 1:ch:nC*Ko
    cols = c0:min(nC*Ko,c0+ch-1);
    CE(:,cols) = A*(b36_gp_k(gp,Zev,C.Z(cols,:)) - Wev'*Wc(:,cols));
end
% candidate observation blocks (latent covariance + nugget)
Sig = zeros(Ko,Ko,nC);
for k1 = 1:Ko
    r1 = (0:nC-1)'*Ko + k1;
    for k2 = k1:Ko
        r2 = (0:nC-1)'*Ko + k2;
        ck = b36_gp_kpair(gp,C.Z(r1,:),C.Z(r2,:)) - sum(Wc(:,r1).*Wc(:,r2),1)';
        Sig(k1,k2,:) = reshape(ck,1,1,nC); Sig(k2,k1,:) = reshape(ck,1,1,nC);
    end
end
Sig = Sig + tau2*eye(Ko);
s2E = zeros(nE,1); Ucand = zeros(nC*Ko,0); valid = true(nC,1);
for b = 1:B
    base = b36_phi2gain(al,min(s2E./v,1-1e-12));
    gain = -Inf(nC,1);
    for cc = find(valid)'
        cols = (cc-1)*Ko + (1:Ko);
        Lc = cholj(Sig(:,:,cc));
        U = CE(:,cols)/Lc';
        rho = min((s2E + sum(U.^2,2))./v,1-1e-12);
        gain(cc) = EW'*(b36_phi2gain(al,rho) - base);
    end
    [score(b),cs] = max(gain);
    if ~isfinite(score(b)), sel = sel(1:b-1); break; end
    sel(b) = cs; valid(cs) = false;
    % ---- exact Schur update on the selected block
    cols = (cs-1)*Ko + (1:Ko);
    Kc = b36_gp_k(gp,C.Z,C.Z(cols,:)) - Wc'*Wc(:,cols);
    if ~isempty(Ucand), Kc = Kc - Ucand*Ucand(cols,:)'; end
    Ls = cholj(Kc(cols,:) + tau2*eye(Ko));
    Ue = CE(:,cols)/Ls'; Un = Kc/Ls';
    CE = CE - Ue*Un';
    s2E = s2E + sum(Ue.^2,2);
    Ucand = [Ucand, Un]; %#ok<AGROW>
    for cc = find(valid)'
        r = (cc-1)*Ko + (1:Ko);
        Sig(:,:,cc) = Sig(:,:,cc) - Un(r,:)*Un(r,:)';
    end
end
info.scores = score(1:numel(sel)); info.explainedEventFraction = sum(EW.*min(s2E./v,1))/sum(EW);
info.s2E = s2E; info.v = v; info.EW = EW;
end

function L = cholj(S)
S = (S+S')/2; j = 0; base = max(mean(diag(S)),realmin);
for it = 1:8
    [L,p] = chol(S + j*eye(size(S)),'lower');
    if p==0, return; end
    if j==0, j = 1e-12*base; else, j = 10*j; end
end
error('b36_acq_sur:chol','block covariance not positive definite');
end
