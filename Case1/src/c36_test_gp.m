function c36_test_gp()
%C36_TEST_GP  Fixed-data accuracy and calibration of the surrogate variants.
c = c36_config(); p = c.pGrid;
bank = c36_bank(4096,77,c.tGrid);
Eref = c36_true_envelope(bank.xi,p,c.tGrid);
dstar = [0.5439 1.1137];
pfTrue = c36_pf(Eref,p,dstar,bank.zeta);
h0 = c36_margin(Eref,p,dstar,bank.zeta);
crit = abs(h0) < 0.15;
specs = {struct('spatial','op','mean','lf1'), struct('spatial','op','mean','lf1','nonstat',true), ...
         struct('spatial','op','mean','const'), struct('spatial','op','mean','const','nonstat',true), ...
         struct('spatial','se','mean','const')};
names = {'OP+LFres','OP+LFres+NS','OP only','OP+NS','SE (OK)'};
% held-out states: population states (realistic) and high-load tail states
rs = RandStream('mt19937ar','Seed',5);
St = bank.S(rs.randperm(size(bank.S,1),400),:);
Stail = St(sqrt(St(:,1).^2+St(:,2).^2) > 2,:);
Yt = c36_loadeffect(St,p,'hf'); Ytl = c36_loadeffect(Stail,p,'hf');
fprintf('%d tail states of 400\n',size(Stail,1));
for n = [20 40 80]
    S = c36_lhs_states(n,11); Y = c36_loadeffect(S,p,'hf');
    for k = 1:numel(specs)
        tic; m = c36_cgp_fit(S,Y,specs{k}); tf = toc;
        Emu = envelope(m,bank);
        pfS = c36_pf(Emu,p,dstar,bank.zeta);
        h1 = c36_margin(Emu,p,dstar,bank.zeta);
        [mt,vt] = c36_cgp_predict(m,St); zs = (Yt - mt)./sqrt(max(vt,1e-300));
        [ml,vl] = c36_cgp_predict(m,Stail); zl = (Ytl - ml)./sqrt(max(vl,1e-300));
        % calibration of the maximum over p (the quantity entering the margin)
        [~,jm] = max(mt,[],1); ix = sub2ind(size(mt),jm,1:size(mt,2));
        zmax = (Yt(ix) - mt(ix))./sqrt(vt(ix));
        [~,jl] = max(ml,[],1); il = sub2ind(size(ml),jl,1:size(ml,2));
        zlmax = (Ytl(il) - ml(il))./sqrt(vl(il));
        g = ''; if strcmp(specs{k}.spatial,'op') && isfield(m.comp{1},'gamma'), g = sprintf(' gamma=[%s]',num2str(cellfun(@(x) x.gamma,m.comp),3)); end
        fprintf('N=%3d %-12s Pf rel %7.1f%% hRMSE=%.4f crit=%.4f | z-sd all %.2f tail %.2f | peak z-sd all %.2f tail %.2f cover95 tail %.2f | fit %.2fs%s\n', ...
            n,names{k},100*(pfS-pfTrue)/pfTrue,rms(h1-h0),rms(h1(crit)-h0(crit)),std(zs(:)),std(zl(:)), ...
            std(zmax),std(zlmax),mean(abs(zlmax)<1.96),tf,g);
    end
end
end

function E = envelope(m,bank)
np = numel(m.pGrid); N = bank.N; nt = bank.nt; E = -inf(N,np);
for it = 1:nt
    rows = (it-1)*N + (1:N);
    mu = c36_cgp_predict(m,bank.S(rows,:));
    E = max(E,mu');
end
end
