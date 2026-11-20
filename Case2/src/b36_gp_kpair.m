function k = b36_gp_kpair(gp,A,B)
%B36_GP_KPAIR  Prior covariance of matching rows, k(A_i,B_i) (column vector).
R2 = zeros(size(A,1),1);
for g = 1:numel(gp.ell)
    R2 = R2 + sum(((A(:,gp.grp==g) - B(:,gp.grp==g))/gp.ell(g)).^2,2);
end
r = sqrt(5*R2); k = gp.sf2*(1 + r + r.^2/3).*exp(-r);
if gp.gamma > 0
    k = k.*sqrt(1 + gp.gamma*A(:,end).^2).*sqrt(1 + gp.gamma*B(:,end).^2);
end
end
