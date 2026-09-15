function summary = aggregateGridSearchResults()
% aggregateGridSearchResults.m
% Scans results/gridSearch/task_*.mat (each written the moment its task
% finishes, by Run_GridSearch.m's parfor loop) and rebuilds the summary
% table -- safe to run at ANY time, including while Run_GridSearch.m is
% still running elsewhere, or after it was interrupted partway through.
% Prints how many of the expected tasks are present vs. still missing.
% Writes results/gridSearch/summary.csv (human-readable) and
% results/gridSearch/summary.mat (the same table).
%
% Calls: none
% Called by: Run_GridSearch.m (automatically, after its parfor loop);
%   also fine to call directly/standalone at any time to check progress.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','gridSearch');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);

if isempty(files)
    summary = table();
    warning('aggregateGridSearchResults:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    if isfield(o,'astarIterMax')
        astarIterMax = o.astarIterMax; %#ok<*PROPLC>
    else
        astarIterMax = NaN; % saved before this field was added (e.g. the 16 EHA tasks from the first pass)
    end
    rows{i} = table({o.label}, o.seaStateIdx, o.Hs, o.Tp, o.probability, ...
        o.bestPressure/1e6, o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, astarIterMax, ...
        'VariableNames', {'Drivetrain','SeaStateIdx','Hs_m','Tp_s','Probability', ...
        'BestPressure_MPa','MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW','AstarIterMax'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'Drivetrain','SeaStateIdx'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');

fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

if any(~isnan(summary.AstarIterMax))
    fprintf(['\nNOTE: DHD rows (AstarIterMax column) were run with the A* search''s\n' ...
        'per-control-window iteration budget capped (see control/MPC_Astar.m and\n' ...
        'Run_GridSearch.m''s comments) to make this sweep tractable -- a single DHD4\n' ...
        'simulation was taking many hours uncapped. These numbers reflect that\n' ...
        'tighter search budget, not the algorithm''s full potential.\n']);
end

% Report what's still missing, if this is being run mid-sweep
expectedLabels = {'PassivePump','DHD2','DHD3','DHD4','EHA_mech','EHA_elec'};
nSeaStates = max(summary.SeaStateIdx);
fprintf('\n--- Completion check (assuming %d sea states) ---\n', nSeaStates);
for L = expectedLabels
    have = sum(strcmp(summary.Drivetrain, L{1}));
    fprintf('%-12s %d/%d\n', L{1}, have, nSeaStates);
end
end
