function fe = c2_assemble(d,xi,cfg)
%C2_ASSEMBLE Hermite-cubic beam elements: consistent mass and the full
%   variable-EI weak form (no response reduction).
%   loadMap*qGauss maps the load density to consistent nodal forces. Only the
%   end displacements are constrained (pinned ends); zero end moment is the
%   natural boundary condition of the weak form, so the rotations are free.
if nargin<3, cfg=c2_config(); end
ne=cfg.nElem; nd=2*(ne+1); le=cfg.L/ne;
ng=cfg.nGauss;
assert(ng>=3 && ng==round(ng),'c2:Quadrature','nGauss must be an integer >= 3');
% Golub-Welsch: Gauss-Legendre quadrature rule (not a reduced-order model).
j=(1:ng-1)'; beta=j./sqrt(4*j.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[z,order]=sort(diag(D)); wt=2*V(1,order)'.^2;
r=(z+1)/2; wt=wt*le/2; [N,~,B]=c2_hermite(r,le);
qX=reshape(((0:ne-1)*le+r*le),[],1);
sec=c2_section(d,xi,qX,cfg);
ii=zeros(ne*16,1); jj=ii; mm=ii; kk=ii;
li=zeros(ne*4*ng,1); lj=li; lv=li; qw=repmat(wt,ne,1);
for e=1:ne
    dof=(2*e-1):(2*e+2); iq=(e-1)*ng+(1:ng);
    Me=N'*(N.*(wt.*sec.rhoA(iq)));
    Ke=B'*(B.*(wt.*sec.EI(iq)));
    [ir,ic]=ndgrid(dof,dof); ix=(e-1)*16+(1:16);
    ii(ix)=ir(:); jj(ix)=ic(:); mm(ix)=Me(:); kk(ix)=Ke(:);
    [ir,ic]=ndgrid(dof,iq); ix=(e-1)*4*ng+(1:4*ng);
    values=N'.*wt'; li(ix)=ir(:); lj(ix)=ic(:); lv(ix)=values(:);
end
M=sparse(ii,jj,mm,nd,nd); K=sparse(ii,jj,kk,nd,nd);
M=(M+M')/2; K=(K+K')/2;
fixed=[1,2*ne+1]; free=setdiff(1:nd,fixed);
fe=struct('M',M,'K',K,'C',cfg.dampingGamma*M,'free',free, ...
    'fixed',fixed,'nDof',nd,'nElem',ne,'le',le,'qX',qX, ...
    'qW',qw,'loadMap',sparse(li,lj,lv,nd,ne*ng));
end
