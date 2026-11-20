function out = b36_hf(d,xi,mode)
%B36_HF  One high-fidelity history: 80 Hermite-cubic elements (c2_assemble),
%   C = 0.8 M, average-acceleration Newmark with dt = 1 ms on 0-6 s, written
%   exactly in the full modal basis (no truncation), plus the retained working block.
%   out.W    : 1-by-nW window maxima of |sigma| on the working grid  [MPa]
%              (the observation block; trajectory working peak = max(out.W))
%   out.peak : max_{x,t}|sigma| on the full HF output grid (400 x 6001) [MPa]
%              (only for mode 'full'; used by the independent validation)
%   mode 'work' (default) : working block only (same FE solution, the stress
%                           recovery is restricted to the retained grid)
%        'full'           : also the full-grid peak
%        'field'          : also the working field out.Sw
%   One call = one HF unit for every method.
if nargin<3 || isempty(mode), mode = 'work'; end
if islogical(mode), if mode, mode = 'field'; else, mode = 'work'; end, end
c = b36_config(); cfg = c.hf;
fe = c2_assemble(d,xi,cfg); f = fe.free;
M = full(fe.M(f,f)); K = full(fe.K(f,f));
[V,lam] = eig(K,M,'vector'); [lam,id] = sort(real(lam)); V = real(V(:,id));
V = V./sqrt(sum(V.*(M*V),1));
t = 0:cfg.dt:cfg.T; F = zeros(numel(f),numel(t));
for k = 1:cfg.loadChunkSize:numel(t)
    ix = k:min(k+cfg.loadChunkSize-1,numel(t));
    F(:,ix) = fe.loadMap(f,:)*c2_load(fe.qX,t(ix),xi,cfg);
end
r = V'*F; q = zeros(size(r)); h = cfg.dt; g = cfg.dampingGamma;
for j = 1:numel(lam)
    den = [4+2*g*h+lam(j)*h^2, -8+2*lam(j)*h^2, 4-2*g*h+lam(j)*h^2];
    q(j,:) = filter(h^2*[1 2 1],den,r(j,:));
end
opw = c2_field_operators(d,xi,cfg,c.xw);
Sw = abs((opw.S(:,f)*V)*q(:,c.twIdx))/1e6;
W = accumarray(c.winOf(:),max(Sw,[],1)',[c.nW 1],@max)';
out = struct('W',W,'peakW',max(W),'omega',sqrt(lam(1:4))','peak',NaN,'tPeak',NaN,'xPeak',NaN);
if any(strcmp(mode,{'full','field'}))
    op = c2_field_operators(d,xi,cfg);
    sig = (op.S(:,f)*V)*q;
    [pk,ix] = max(abs(sig),[],1); [gp,it] = max(pk);
    out.peak = gp/1e6; out.tPeak = t(it); out.xPeak = op.x(ix(it)); out.signPeak = sign(sig(ix(it),it));
end
if strcmp(mode,'field'), out.Sw = Sw; end
end
