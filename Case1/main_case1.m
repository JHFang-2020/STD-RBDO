%% MAIN_CASE1  One-click script for Case 1 (analytical benchmark): Tables 3 and 4 of the paper.
%
%   Choose MODE and NWORKERS below and run the script (F5, or "run main_case1" in this folder).
%
%     MODE = 'test'    unit tests of the surrogate and of the event covariances (about 1 min)
%     MODE = 'quick'   run 1 of every configuration of Tables 3 and 4 (14 runs), then both tables
%     MODE = 'full'    10 runs of every configuration, as in the paper (7 x 10 runs for Table 3,
%                      7 further configurations x 10 runs for Table 4), then both tables
%     MODE = 'tables'  rebuild Tables 3 and 4 from the runs already in results/main (no new runs;
%                      the first call builds the exact-model scoring data in results/truth, about 7 min)
%
%   Outputs (folder results/, created on first use):
%     results/main/*.mat           one file per run (method, surrogate, run number)
%     results/reference_rows.mat   first three rows of Table 3
%     results/table3_case1.csv     Table 3, optimized designs of Case 1
%     results/table4_case1.csv     Table 4, ablation study of Case 1
%     results/truth/               cached exact-model data (used only to score designs)
%
%   Finished runs are skipped, so the script can be stopped and restarted at any time.
%   Requirements: MATLAB, Optimization Toolbox, Statistics and Machine Learning Toolbox;
%   Parallel Computing Toolbox optional (without it the runs are executed one after another).

MODE     = 'full';   % 'test' | 'quick' | 'full' | 'tables'
NWORKERS = 4;        % number of parallel workers

%% ------------------------------------------------------------------------------------------
case1Dir = fileparts(mfilename('fullpath'));
addpath(case1Dir, fullfile(case1Dir,'src'));
case1Timer = tic;
switch lower(MODE)
    case 'test'
        c36_test_gp();
        c36_test_events();
    case 'quick'
        run_case1_table3(1, NWORKERS);
        run_case1_table4(1, NWORKERS);
    case 'full'
        run_case1_table3(1:10, NWORKERS);
        run_case1_table4(1:10, NWORKERS);
    case 'tables'
        runFiles = dir(fullfile(case1Dir,'results','main','*_s*.mat'));
        runTok = regexp({runFiles.name},'_s(\d+)(_[A-Za-z0-9]+)?\.mat$','tokens','once');
        runSeeds = unique(cellfun(@(t) str2double(t{1}), runTok(~cellfun(@isempty,runTok))));
        if isempty(runSeeds)
            error('main_case1:noRuns','No runs in %s. Use MODE = ''quick'' or ''full'' first.', ...
                fullfile(case1Dir,'results','main'));
        end
        if ~isfile(fullfile(case1Dir,'results','reference_rows.mat')), case1_reference_rows(); end
        make_case1_table3(runSeeds);
        make_case1_table4(runSeeds);
    otherwise
        error('main_case1:mode','Unknown MODE ''%s''. Use ''test'', ''quick'', ''full'' or ''tables''.',MODE);
end
fprintf('\nmain_case1 (MODE = %s) finished in %.1f min.\n', MODE, toc(case1Timer)/60);
