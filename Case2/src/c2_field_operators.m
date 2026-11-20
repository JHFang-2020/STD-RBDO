function op = c2_field_operators(d,xi,cfg,x,element)
%C2_FIELD_OPERATORS Linear maps from the global displacement DOFs to the
%   deflection w, the rotation and the signed bending stress.
%   Default points: equally spaced in every element, with the duplicated
%   left/right nodes kept. If x is given but element is not, an interior node
%   is evaluated in the element to its right and x = L in the last element.
%   Curvature peaks are therefore not averaged or smoothed.
ne=cfg.nElem; le=cfg.L/ne; nd=2*(ne+1);
if nargin<4 || isempty(x)
    rr=linspace(0,1,cfg.samplePtsPerElem)';
    x=reshape((0:ne-1)*le+rr*le,[],1);
    element=repelem((1:ne)',numel(rr));
else
    x=x(:); assert(all(x>=0 & x<=cfg.L));
    if nargin<5 || isempty(element)
        element=min(floor(x/le+1e-10)+1,ne);
    end
end
element=element(:); r=(x-(element-1)*le)/le;
[N,Nx,Nxx]=c2_hermite(r,le); sec=c2_section(d,xi,x,cfg);
rows=repelem((1:numel(x))',4); cols=reshape((2*element-1+(0:3))',[],1);
W=sparse(rows,cols,reshape(N',[],1),numel(x),nd);
R=sparse(rows,cols,reshape(Nx',[],1),numel(x),nd);
S=sparse(rows,cols,reshape((Nxx.*(-0.5*sec.E.*sec.h))',[],1),numel(x),nd);
op=struct('x',x,'element',element,'W',W,'R',R,'S',S);
end
