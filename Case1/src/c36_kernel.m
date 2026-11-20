function K = c36_kernel(X1,X2,ell,type)
%C36_KERNEL  Unit-variance stationary kernels with ARD length scales.
%   type 'm52' : Matern-5/2 (load-state kernel), 'se' : squared exponential.
X1 = X1./ell(:)'; X2 = X2./ell(:)';
r2 = max(sum(X1.^2,2) + sum(X2.^2,2)' - 2*(X1*X2'),0);
switch type
    case 'm52'
        r = sqrt(5*r2);
        K = (1 + r + r.^2/3).*exp(-r);
    case 'se'
        K = exp(-0.5*r2);
    otherwise
        error('c36_kernel:type','unknown kernel');
end
end
