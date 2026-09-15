function summary = aggregateAstarHorizonResults()
% aggregateAstarHorizonResults.m
% Scans results/astarHorizon/task_*.mat and rebuilds the summary table.
% Safe to run at any time, including mid-study.
% Calls: none
% Called by: Run_AstarHorizonSearch.m (automatically); fine standalone too.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','astarHorizon');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed A*-horizon task file(s) in %s\n', numel(files), resultsDir);

if isempty(files)
    summary = table();
    warning('aggregateAstarHorizonResults:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table(o.seaStateIdx, o.Hs, o.Tp, o.mAstar, o.highPressure/1e6, ...
        o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, ...
        'VariableNames', {'SeaStateIdx','Hs_m','Tp_s','mAstar','Pressure_MPa', ...
        'MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'SeaStateIdx','mAstar'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s\n', fullfile(resultsDir,'summary.csv'));
disp(summary);
end
