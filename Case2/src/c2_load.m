function [q,history] = c2_load(x,t,xi,cfg)
%C2_LOAD Two moving loads, normalised on the finite span; output on the
%   Cartesian grid numel(x)-by-numel(t).
%   The whole history uses the same 12 load coefficients xi(6:17); the function
%   neither reads nor changes the random-number state.
%   Order within each group: [cos1 sin1 cos2 sin2 cos3 sin3].
if nargin<4, cfg=c2_config(); end
validateattributes(xi,{'numeric'},{'vector','numel',17,'finite'});
x=x(:); t=t(:)'; xi=xi(:);
assert(all(x>=0 & x<=cfg.L) && all(t>=0 & t<=cfg.T));
a=[1 0.6 0.35]/sqrt(1+0.6^2+0.35^2);
phase=(2*pi*(1:3)'/cfg.TL)*t;
C=cos(phase); S=sin(phase); Y=zeros(2,numel(t));
for r=1:2
    co=reshape(xi(6+6*(r-1):11+6*(r-1)),2,3);
    Y(r,:)=(a.*co(1,:))*C+(a.*co(2,:))*S;
end
Z=[Y(1,:);cfg.loadCorrelation*Y(1,:)+sqrt(1-cfg.loadCorrelation^2)*Y(2,:)];
c=cfg.L*[0.5+0.35*tanh(1.1*sin(2*pi*t/cfg.TL)+0.30*Z(1,:)); ...
    0.5-0.35*tanh(1.1*sin(2*pi*t/cfg.TL+0.7)+0.30*Z(2,:))];
P=zeros(2,numel(t)); envelope=P; q=zeros(numel(x),numel(t));
for r=1:2
    envelope(r,:)=sin(pi*t/cfg.TL).^2.* ...
        (1+(-1)^(r+1)*0.35*sin(2*pi*t/cfg.TL));
    envelope(r,t>=cfg.TL)=0; envelope(r,t==0)=0;
    P(r,:)=cfg.P0(r)*envelope(r,:).* ...
        exp(cfg.loadLogScale*Z(r,:)-0.5*cfg.loadLogScale^2).* ...
        (1+0.45*sin(cfg.omega(r)*t+cfg.phi(r)));
    ell=cfg.ell(r);
    denom=ell*sqrt(pi/2)*(erf((cfg.L-c(r,:))/(sqrt(2)*ell))+ ...
        erf(c(r,:)/(sqrt(2)*ell)));
    q=q+exp(-(x-c(r,:)).^2/(2*ell^2)).*(P(r,:)./denom);
end
history=struct('t',t,'Y',Y,'Z',Z,'centers',c,'P',P,'envelope',envelope);
end
