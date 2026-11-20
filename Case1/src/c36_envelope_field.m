function E = c36_envelope_field(F)
%C36_ENVELOPE_FIELD  Time envelope max_t F(p,n,t) with 3-point parabolic refinement.
%   F : np-by-N-by-nt field on a uniform time grid;  E : N-by-np
[np,N,nt] = size(F);
[fmax,j] = max(F,[],3);                 % np-by-N
E = fmax;
in = j>1 & j<nt;
if any(in(:))
    [ip,in2] = find(in);
    jj = j(in);
    lin0 = sub2ind([np N nt],ip,in2,jj-1);
    lin2 = sub2ind([np N nt],ip,in2,jj+1);
    f0 = double(F(lin0)); f1 = double(fmax(in)); f2 = double(F(lin2));
    den = f0 - 2*f1 + f2;              % <= 0 at a maximum
    ok = den < -1e-14;
    delta = zeros(size(f1));
    delta(ok) = 0.5*(f0(ok)-f2(ok))./den(ok);
    delta = max(min(delta,0.5),-0.5);
    E(in) = f1 - 0.25*(f0-f2).*delta;
end
E = double(E');
end
