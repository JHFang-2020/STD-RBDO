function L = b36_lf_only(variants)
%B36_LF_ONLY  Designs obtained from the simplified physics alone (no HF, no
%   learning), scored with the exact working-population truth: shows what
%   each physics level gives without the discrepancy GP.
%   Same structure as the exact optimum (d1 at its lower, d2 at its upper
%   bound; cheapest d3 by a monotone root of the LF 41st-largest peak).
%   variants : cell array of simplified-physics variants (default
%   {'Efield','nomE','qs'}); Table 5 uses 'Efield'.
if nargin<1 || isempty(variants), variants = {'Efield','nomE','qs'}; end
c = b36_config(); P = b36_pop('work'); TR = b36_truth_ref();
[P.cache,~] = b36_cache(P.cache,P.N,P.xi,(1:P.N)');
kst = floor(c.target*P.N) + 1; rows = {};
for v = variants(:)'
    lfo = struct('variant',v{1});
    g = @(d3) qLF([c.lower(1) c.upper(2) d3]) - c.allow;
    if g(c.lower(3)) <= 0, d3 = c.lower(3);
    elseif g(c.upper(3)) > 0, d3 = c.upper(3);
    else, d3 = fzero(g,[c.lower(3) c.upper(3)],optimset('TolX',1e-6));
    end
    d = [c.lower(1) c.upper(2) d3];
    [P,S] = b36_screen(P,d,0.3,struct());                 % truth screening with the reference LF
    T = b36_truth_at(P,S,struct('theta',0.25,'maxHF',3000));
    rows{end+1} = table(string(v{1}),d,sum(d),(sum(d)-TR.CW)/TR.CW,T.q,T.qErr,T.nFail,T.pf,T.exact, ...
        'VariableNames',{'lf','d','cost','costGapW','qTrue','qErr','nFail','pfTrueW','exact'}); %#ok<AGROW>
end
L = vertcat(rows{:});
writetable(L,fullfile(fileparts(mfilename('fullpath')),'..','results','C2_lf_only_designs.csv'));
disp(L);
    function q = qLF(d)
        [~,ld] = b36_cache(P.cache,P.N,P.xi,(1:P.N)');
        o = b36_lf(d,P.xi,ld,lfo); pk = sort(max(o.W,[],2),'descend'); q = pk(kst);
    end
end
