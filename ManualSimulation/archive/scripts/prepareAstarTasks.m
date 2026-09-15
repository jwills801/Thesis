function prepareAstarTasks()
% prepareAstarTasks.m
% Builds the A*-horizon-length study's task list and writes it to
% results/astarHorizon/taskList.csv, one row per task. Deliberately
% scoped down (see run_overnight.sh's header for the full rationale):
%   - DHD, 2 rails only
%   - 3 representative sea states (indices 1, 5, 8 of the 8 Humboldt bins)
%   - horizon lengths {2,3,4,5,6,7,8} (default elsewhere is 5)
%   - pressure held FIXED per sea state at whatever
%     results/gridSearch/summary.csv already found best for DHD2 there
% Each row is run as its OWN OS process by run_overnight.sh
% (runSingleAstarTask.m does the actual work for one row) so a slow
% horizon length can be killed with a real timeout -- this machine has no
% Parallel Computing Toolbox, so there's no in-MATLAB way to cancel a
% running computation.
%
% Calls: wave/humboldtSeaStates.m
% Called by: run_overnight.sh (via matlab -batch)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'wave'));

resultsDir = fullfile(repoRoot,'results','astarHorizon');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

seaStates = humboldtSeaStates();
seaStateIdxs = [1, 5, 8];
horizonLengths = [2 3 4 5 6 7 8];

summaryFile = fullfile(repoRoot,'results','gridSearch','summary.csv');
if ~isfile(summaryFile)
    error('prepareAstarTasks:missingSummary','results/gridSearch/summary.csv not found -- run Run_GridSearch.m first.');
end
summary = readtable(summaryFile);

TaskIdx = []; SeaStateIdx = []; Hs = []; Tp = []; HighPressure_Pa = []; MAstar = [];
for s = seaStateIdxs
    rowMask = strcmp(summary.Drivetrain,'DHD2') & summary.SeaStateIdx==s;
    if ~any(rowMask)
        warning('prepareAstarTasks:noPressureResult','No DHD2 result for sea state %d, skipping.',s);
        continue
    end
    bestPressureMPa = summary.BestPressure_MPa(find(rowMask,1));
    for m = horizonLengths
        TaskIdx(end+1,1) = numel(TaskIdx)+1; %#ok<AGROW>
        SeaStateIdx(end+1,1) = s; %#ok<AGROW>
        Hs(end+1,1) = seaStates.Hs(s); %#ok<AGROW>
        Tp(end+1,1) = seaStates.Tp(s); %#ok<AGROW>
        HighPressure_Pa(end+1,1) = bestPressureMPa*1e6; %#ok<AGROW>
        MAstar(end+1,1) = m; %#ok<AGROW>
    end
end

taskList = table(TaskIdx, SeaStateIdx, Hs, Tp, HighPressure_Pa, MAstar);
writetable(taskList, fullfile(resultsDir,'taskList.csv'));
fprintf('Wrote %d A*-horizon tasks to %s\n', height(taskList), fullfile(resultsDir,'taskList.csv'));
end
