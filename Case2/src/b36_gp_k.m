function K = b36_gp_k(gp,A,B)
%B36_GP_K  Prior covariance sf2 * w(A) w(B)' .* k_M52(A,B) of the discrepancy GP.
R2 = zeros(size(A,1),size(B,1));
for g = 1:numel(gp.ell)
    a = A(:,gp.grp==g)/gp.ell(g); b = B(:,gp.grp==g)/gp.ell(g);
    R2 = R2 + max(sum(a.^2,2) + sum(b.^2,2)' - 2*(a*b'),0);
end
r = sqrt(5*R2); K = gp.sf2*(1 + r + r.^2/3).*exp(-r);
if gp.gamma > 0
    K = K.*(sqrt(1 + gp.gamma*A(:,end).^2)*sqrt(1 + gp.gamma*B(:,end).^2)');
end
end
