function D = c36_loadeffect(S,p,model)
%C36_LOADEFFECT  Load-effect field D(s,p) = A - 1.8(p-c)^2 - 4(p-c)^4.
%   S : M-by-3 states, p : column of spatial nodes (np-by-1)
%   D : np-by-M (one column = one spatial snapshot)
%   model 'hf' is the expensive simulator (one column = one HF call);
%   model 'lf' is the simplified physics used as GP prior mean.
if nargin<3, model = 'hf'; end
c = c36_config();
[A,cen] = c36_fields(S,model);
r = p(:) - cen(:)';               % np-by-M
D = A(:)' - c.k*r.^2 - c.lambda*r.^4;
end
