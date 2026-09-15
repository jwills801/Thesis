function summary = aggregateFinalGridSearchResults()
% aggregateFinalGridSearchResults.m
% Scans results/gridSearchFinal/task_*.mat (each written the moment its
% task finishes) and rebuilds summary.csv/.mat -- safe to run at any
% time, including mid-sweep or after an interrupted session.
%
% Calls: none
% Called by: Run_FinalGridSearch_100ms.m (automatically); also fine to
%   call directly/standalone at any time to check progress.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','gridSearchFinal');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);
if isempty(files)
    summary = table();
    warning('aggregateFinalGridSearchResults:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    mAstarVal = o.mAstar; if isnan(mAstarVal), mAstarVal = NaN; end
    rows{i} = table({o.label}, o.seaStateIdx, o.Hs, o.Tp, o.probability, ...
        o.bestPressure/1e6, o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, ...
        o.nAstarCapHits, mAstarVal, ...
        'VariableNames', {'Drivetrain','SeaStateIdx','Hs_m','Tp_s','Probability', ...
        'BestPressure_MPa','MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW', ...
        'AstarCapHits','MAstar'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'Drivetrain','SeaStateIdx'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

if any(summary.AstarCapHits > 0)
    fprintf('\nWARNING: some DHD tasks hit the A* iteration cap (astarIterMax) -- check which and investigate.\n');
end

expectedLabels = {'PassivePump','DHD2','DHD3','DHD4','EHA_mech','EHA_elec'};
nSeaStates = max(summary.SeaStateIdx);
fprintf('\n--- Completion check (assuming %d sea states) ---\n', nSeaStates);
for L = expectedLabels
    have = sum(strcmp(summary.Drivetrain, L{1}));
    fprintf('%-12s %d/%d\n', L{1}, have, nSeaStates);
end
end
