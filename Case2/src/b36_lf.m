function out = b36_lf(d,xi,ld,opts)
%B36_LF  Simplified physics of the histories xi (rows) at design d.
%   sigma = 6 M / (b h(x)^2)                 (exact section relation)
%   M     = Mqs + sum_j Psi_j(x) eta_j(t)    (mode acceleration)
%   Mqs   : exact quasi-static moment (b36_loads, operator level)
%   eta_j : dynamic excess q_j - f_j/k_j of J assumed sine modes of the
%           history's own section and modulus field (single-term
%           Rayleigh-Ritz per mode, no modal coupling), Newmark average
%           acceleration on the 5 ms working step.
%   opts.variant 'Efield' (default) | 'nomE' (nominal modulus) | 'qs' (no dynamics)
%   out.W   : n-by-nW window maxima of |sigma_LF| [MPa]
%   out.Wqs : n-by-nW window maxima of the quasi-static stress (exact) [MPa]
%   out.S   : nx-by-ntw-by-n field (only if opts.field)
if nargin<4, opts = struct(); end
n = size(xi,1); ch = 64;
if n > ch && ~getf(opts,'field',false)
    % bounded memory: process blocks of histories, keep only window maxima
    out = struct('W',zeros(n,0),'Wqs',zeros(n,0),'h',zeros(n,0));
    for k0 = 1:ch:n
        ii = k0:min(n,k0+ch-1);
        ldk = struct('Mqs',ld.Mqs(:,:,ii),'f',ld.f(ii,:,:),'n',numel(ii));
        o = b36_lf(d,xi(ii,:),ldk,opts);
        out.W(ii,1:size(o.W,2)) = o.W; out.Wqs(ii,1:size(o.Wqs,2)) = o.Wqs; out.h(ii,1:size(o.h,2)) = o.h;
    end
    return;
end
variant = getf(opts,'variant','Efield');
c = b36_config(); cfg = c.hf; L = c.L; b = c.b; J = c.J;
xq = c.xq; wq = trapzw(xq);
X = d(:)'.*exp(c.sigmaH*xi(:,1:3) - 0.5*c.sigmaH^2);       % n-by-3
H = X*bern(c.xw/L)';  Hq = X*bern(xq/L)';                % n-by-nx, n-by-nq
Mqs = double(ld.Mqs);
den = reshape((b*H.^2)',c.nx,1,n);
Sqs = abs(6*Mqs./den)/1e6;
if strcmp(variant,'qs')
    S = Sqs;
else
    k = pi/L;
    if strcmp(variant,'Efield')
        co = cos(k*xq');
        Eq = cfg.E0*exp(cfg.sigmaE*(sqrt(0.7)*xi(:,4) + sqrt(0.3)*xi(:,5).*co) - 0.5*cfg.sigmaE^2*(0.7+0.3*co.^2));
        cox = cos(k*c.xw');
        Ex = cfg.E0*exp(cfg.sigmaE*(sqrt(0.7)*xi(:,4) + sqrt(0.3)*xi(:,5).*cox) - 0.5*cfg.sigmaE^2*(0.7+0.3*cox.^2));
    else
        Eq = cfg.E0*ones(n,numel(xq)); Ex = cfg.E0*ones(n,c.nx);
    end
    Iq = b*Hq.^3/12; Aq = b*Hq; Ix = b*H.^3/12;
    M = Mqs;
    for j = 1:J
        sj = sin(j*pi*xq'/L).^2;
        mj = (cfg.rho*Aq.*sj)*wq;                         % n-by-1 modal mass
        kj = (Eq.*Iq.*sj)*wq*(j*pi/L)^4;                  % modal stiffness
        fj = double(ld.f(:,:,j));
        qj = newmark(fj./mj,kj./mj,c.gammaD,c.dtw);
        eta = qj(:,c.twInLF) - fj(:,c.twInLF)./kj;        % dynamic excess
        Psi = Ex.*Ix*(j*pi/L)^2.*sin(j*pi*c.xw'/L);       % sagging moment of w = sin (w along the load)
        M = M + reshape(Psi',c.nx,1,n).*reshape(eta',1,c.ntw,n);
    end
    S = abs(6*M./den)/1e6;
end
out = struct('W',winmax(S,c),'Wqs',winmax(Sqs,c),'h',H);
if getf(opts,'field',false), out.S = S; end
end

function W = winmax(S,c)
% n-by-nW maxima over space and over each 0.1 s time window
smax = squeeze(max(S,[],1));                % ntw-by-n (or ntw-by-1)
if isvector(smax), smax = smax(:); end
n = size(smax,2); W = zeros(n,c.nW);
for w = 1:c.nW
    W(:,w) = max(smax(c.winOf==w,:),[],1)';
end
end

function q = newmark(r,w2,g,h)
% average-acceleration Newmark, zero initial state, vectorised over rows
a1 = 4+2*g*h+w2*h^2; a2 = -8+2*w2*h^2; a3 = 4-2*g*h+w2*h^2;
[n,nt] = size(r); q = zeros(n,nt);
rp = zeros(n,1); rpp = rp; qp = rp; qpp = rp;
for k = 1:nt
    qk = (h^2*(r(:,k) + 2*rp + rpp) - a2.*qp - a3.*qpp)./a1;
    q(:,k) = qk; qpp = qp; qp = qk; rpp = rp; rp = r(:,k);
end
end

function B = bern(s)
s = s(:); B = [(1-s).^2, 2*s.*(1-s), s.^2];
end

function w = trapzw(x)
dx = diff(x(:)); w = zeros(numel(x),1);
w(1:end-1) = w(1:end-1) + dx/2; w(2:end) = w(2:end) + dx/2;
end

function v = getf(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
