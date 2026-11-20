function pop = c36_popeval(model,bank)
%C36_POPEVAL  Posterior-mean load-effect field on all population states.
%   pop.MU : np-by-N-by-nt (single), pop.E : N-by-np time envelope of the mean
c = c36_config(); np = c.np; N = bank.N; nt = bank.nt;
MU = zeros(np,N,nt,'single');
chunkT = max(1,floor(60000/N));
for it0 = 1:chunkT:nt
    its = it0:min(nt,it0+chunkT-1);
    rows = reshape((its-1)*N + (1:N)',[],1);
    mu = c36_cgp_predict(model,bank.S(rows,:));
    MU(:,:,its) = reshape(single(mu),np,N,numel(its));
end
pop = struct('MU',MU,'E',c36_envelope_field(MU),'N',N,'nt',nt);
end
