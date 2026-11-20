function TR = b36_truth_ref(resDir,kind)
%B36_TRUTH_REF  Exact optimum of the working population (HF, never used in learning).
%   kind 'work' (default) or 'val' (independent validation bank, full-grid
%   HF peaks: the reference for the "total error" basis).
%   min sum(d) s.t. q_(k*)(d) < 160 MPa (<=> Pf_W(d) <= target), where q_(k*)
%   is the k*-th largest HF working peak.  q is continuous and decreasing in
%   every d_i, so for given (d1,d2) the cheapest d3 is a 1-D root.  The
%   surrogate studies put the optimum at d1 = lower, d2 = upper bound; this is
%   verified by perturbing both bounds (KKT check) and, if any perturbation is
%   cheaper, by a pattern search in (d1,d2).
if nargin<1 || isempty(resDir), resDir = fullfile(fileparts(mfilename('fullpath')),'..','results','truth'); end
if nargin<2, kind = 'work'; end
if ~isfolder(resDir), mkdir(resDir); end
if strcmp(kind,'val'), f = fullfile(resDir,'truth_ref_val.mat'); else, f = fullfile(resDir,'truth_ref.mat'); end
if isfile(f), s = load(f); TR = s.TR; return; end
c = b36_config();
P = b36_pop(kind); full = strcmp(kind,'val');
if ~full, [P.cache,~] = b36_cache(P.cache,P.N,P.xi,(1:P.N)'); end   % all loads (exact statics)
evlog = {}; nq = 0;
r0 = root([c.lower(1) c.upper(2)]);
best = r0; tab = [c.lower(1) c.upper(2) r0.d3 r0.cost];
for dl = [0.004 0.010 0.020]
    for which = 1:2
        d12 = [c.lower(1) c.upper(2)];
        if which==1, d12(1) = d12(1) + dl; else, d12(2) = d12(2) - dl; end
        r = root(d12); tab = [tab; d12 r.d3 r.cost]; %#ok<AGROW>
        if r.cost < best.cost - 1e-6, best = r; end
    end
end
if best.cost < r0.cost - 1e-6
    % a bound is not active: local pattern search in (d1,d2)
    fo = optimset('Display','off','TolX',1e-4,'TolFun',1e-6,'MaxFunEvals',60);
    x = fminsearch(@(x) costOf(x),best.d12,fo); best = root(min(max(x,c.lower(1:2)),c.upper(1:2)));
end
dW = [best.d12 best.d3];
[P,S] = b36_screen(P,dW,0.3);
T = b36_truth_at(P,S,struct('theta',0.25,'maxHF',6000,'full',full));
if full, pkT = T.peaksFull; else, pkT = T.peaks; end
marg = log(c.allow./pkT);
crit = T.scr(abs(marg) < 0.10); margCrit = marg(abs(marg) < 0.10);
TR = struct('dW',dW,'CW',sum(dW),'massW',c.massCoef*sum(dW),'T',T,'pfW',T.pf,'qW',T.q, ...
    'kstar',T.kstar,'crit',crit,'margCrit',margCrit,'kkt',tab,'nTruthEvals',nq,'log',{evlog});
save(f,'TR');
fprintf('exact W optimum d = (%.5f %.5f %.5f), C = %.5f, Pf_W = %.5f (%d failures), q = %.3f MPa, truth evals %d\n', ...
    dW,sum(dW),T.pf,T.nFail,T.q,nq);

    function v = costOf(x)
        x = min(max(x,c.lower(1:2)),c.upper(1:2)); r = root(x); v = r.cost;
    end
    function r = root(d12)
        g = @(d3) qAt([d12 d3]) - c.allow;
        lo = c.lower(3); hi = c.upper(3);
        if g(hi) > 0, r = struct('d12',d12,'d3',NaN,'cost',Inf); return; end
        if g(lo) <= 0, d3 = lo;
        else, d3 = fzero(g,[lo hi],optimset('TolX',2e-6));
        end
        r = struct('d12',d12,'d3',d3,'cost',sum(d12)+d3);
        fprintf('  root d12=(%.4f %.4f) -> d3=%.5f cost %.5f\n',d12,d3,r.cost);
    end
    function q = qAt(d)
        [P,Sd] = b36_screen(P,d,0.3);
        Td = b36_truth_at(P,Sd,struct('theta',0.12+0.03*full,'full',full,'maxHF',2500));
        if ~Td.exact, Td = b36_truth_at(P,Sd,struct('theta',0.25,'maxHF',6000,'full',full)); end
        q = Td.q; nq = nq + 1;
        evlog{end+1} = [d Td.q Td.nFail Td.nHF Td.maxErr Td.exact];
    end
end
