function [sig,lin,zvAll] = c36_pointsigma(model,bank,P,lin)
%C36_POINTSIGMA  Posterior s.d. of the margin at selected population points.
%   lin : linear indices into the N-by-nt point arrays (= rows of bank.S)
%   sig : s.d. at (s_lin, p_{j(lin)})
if model.lfOnly
    sig = zeros(numel(lin),1); zvAll = zeros(model.r,numel(lin)); return;
end
[~,zvAll] = c36_cw(model,bank.S(lin,:));
Phi2 = model.Phi.^2;                        % np-by-r
sig = sqrt(max(sum(Phi2(P.j(lin),:).*zvAll',2),0));
end
