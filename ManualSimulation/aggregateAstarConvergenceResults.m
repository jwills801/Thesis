function summary = aggregateAstarConvergenceResults()
% aggregateAstarConvergenceResults.m
% Scans results/astarConvergence100ms/task_*.mat (each written the moment
% its task finishes, by Run_AstarConvergenceStudy_100ms.m's parfor loop)
% and rebuilds summary.csv/.mat -- safe to run at ANY time, including
% while the study is still running elsewhere, or after the session
% ended partway through (every completed task's file survives regardless
% of whether this aggregation step itself ever got to run).
%
% Calls: none
% Called by: Run_AstarConvergenceStudy_100ms.m (automatically, at the end
%   of its parfor loop); also fine to call directly/standalone at any time.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','astarConvergence100ms');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);

if isempty(files)
    summary = table();
    warning('aggregateAstarConvergenceResults:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.family}, o.rails, o.mAstar, o.mechRGP, o.elecRGP, ...
        o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.elapsedSec, ...
        'VariableNames', {'Family','Rails','MAstar','MechRGP','ElecRGP', ...
        'AveMechPow_kW','AveElecPow_kW','AstarCapHits','ElapsedSec'});
end
summary = vertcat(rows{:});
summary = sortrows(summary,{'Family','MAstar'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');

fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

% Report what's still missing, if this is being run mid-sweep
expectedFamilies = {'DHD2','DHD3','DHD4'};
fprintf('\n--- Completion check ---\n');
for fm = expectedFamilies
    have = sort(summary.MAstar(strcmp(summary.Family, fm{1})))';
    fprintf('%-6s mAstar values done: %s\n', fm{1}, mat2str(have));
end
end
