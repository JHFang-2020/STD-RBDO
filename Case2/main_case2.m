%% MAIN_CASE2  One-click script for Case 2 (variable-depth beam under moving loads): Table 5 of the paper.
%
%   Choose MODE and NWORKERS below and run the script (F5, or "run main_case2" in this folder).
%
%     MODE = 'test'    unit tests of the surrogate, the margins and the acquisition (about 1 min)
%     MODE = 'quick'   run 1 of every configuration of Table 5 (6 runs), the reference rows, the
%                      validation of the stopped designs, then the table
%     MODE = 'full'    10 runs of every configuration, as in the paper (6 x 10 runs), then the table
%     MODE = 'tables'  rebuild Table 5 from the runs already in results/main (no new runs; missing
%                      reference rows or validations are computed)
%
%   Outputs (folder results/, created on first use):
%     results/main/*.mat               one file per run, e.g. stmsur_EfieldNS_B4_s001.mat
%     results/reference_rows.mat       first three rows of Table 5
%     results/table5_case2.csv         Table 5, optimized designs at stopping in Case 2
%     results/table5_case2_runs.csv    stopped design and validation Pf of every run
%     results/truth/truth_ref.mat      reference solution (exact HF optimum)
%     results/C2_val_cache.mat, C2_ddo.mat, C2_lf_only_designs.csv   cached validations and designs
%
%   Finished runs and validations are skipped, so the script can be stopped and restarted.
%   Requirements: MATLAB, Optimization Toolbox; Parallel Computing Toolbox optional.

MODE     = 'full';   % 'test' | 'quick' | 'full' | 'tables'
NWORKERS = 4;        % number of parallel workers

%% ------------------------------------------------------------------------------------------
case2Dir = fileparts(mfilename('fullpath'));
addpath(case2Dir, fullfile(case2Dir,'src'));
case2Timer = tic;
switch lower(MODE)
    case 'test'
        b36_test_core();
        b36_test_acq();
    case 'quick'
        run_case2_table5(1, NWORKERS);
    case 'full'
        run_case2_table5(1:10, NWORKERS);
    case 'tables'
        runFiles = dir(fullfile(case2Dir,'results','main','*_s*.mat'));
        runTok = regexp({runFiles.name},'_s(\d+)(_[A-Za-z0-9]+)?\.mat$','tokens','once');
        runSeeds = unique(cellfun(@(t) str2double(t{1}), runTok(~cellfun(@isempty,runTok))));
        if isempty(runSeeds)
            error('main_case2:noRuns','No runs in %s. Use MODE = ''quick'' or ''full'' first.', ...
                fullfile(case2Dir,'results','main'));
        end
        make_case2_table5(runSeeds, NWORKERS);
    otherwise
        error('main_case2:mode','Unknown MODE ''%s''. Use ''test'', ''quick'', ''full'' or ''tables''.',MODE);
end
fprintf('\nmain_case2 (MODE = %s) finished in %.1f min.\n', MODE, toc(case2Timer)/60);
