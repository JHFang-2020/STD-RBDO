function sec = c2_section(d,xi,x,cfg)
%C2_SECTION Cross-section and material fields and their analytical
%   derivatives with respect to the physical coordinate x (not s = x/L).
%   Order of xi: manufacturing 1:3, material 4:5, load history Y1 6:11, Y2 12:17.
%   All outputs are numel(x)-by-1; EIx and EIxx enter the full variable-stiffness operator.
if nargin<4, cfg=c2_config(); end
validateattributes(d,{'numeric'},{'vector','numel',3,'finite','positive'});
validateattributes(xi,{'numeric'},{'vector','numel',17,'finite'});
x=x(:); d=d(:); xi=xi(:); s=x/cfg.L;
assert(all(x>=0 & x<=cfg.L),'c2:SpaceRange','x must lie in [0,L]');
X=d.*exp(cfg.sigmaH*xi(1:3)-0.5*cfg.sigmaH^2);
B=[(1-s).^2,2*s.*(1-s),s.^2];
hx=[-2*(1-s),2*(1-2*s),2*s]*X/cfg.L;
hxx=repmat(2*(X(1)-2*X(2)+X(3))/cfg.L^2,size(x));
h=B*X; k=pi/cfg.L; co=cos(k*x); si=sin(k*x);
vE=0.7+0.3*co.^2; vEx=-0.6*k*co.*si;
vExx=-0.6*k^2*(co.^2-si.^2);
VE=sqrt(0.7)*xi(4)+sqrt(0.3)*xi(5)*co;
logE=cfg.sigmaE*VE-0.5*cfg.sigmaE^2*vE;
logEx=-cfg.sigmaE*sqrt(0.3)*xi(5)*k*si-0.5*cfg.sigmaE^2*vEx;
logExx=-cfg.sigmaE*sqrt(0.3)*xi(5)*k^2*co-0.5*cfg.sigmaE^2*vExx;
E=cfg.E0*exp(logE); Ex=E.*logEx; Exx=E.*(logExx+logEx.^2);
A=cfg.b*h; I=cfg.b*h.^3/12;
Ix=cfg.b*h.^2.*hx/4; Ixx=cfg.b*(2*h.*hx.^2+h.^2.*hxx)/4;
sec=struct('x',x,'X',X,'h',h,'hx',hx,'hxx',hxx,'E',E, ...
    'Ex',Ex,'Exx',Exx,'A',A,'I',I,'EI',E.*I, ...
    'EIx',Ex.*I+E.*Ix,'EIxx',Exx.*I+2*Ex.*Ix+E.*Ixx, ...
    'rhoA',cfg.rho*A,'expectedMass',cfg.rho*cfg.b*cfg.L*sum(d)/3);
end
