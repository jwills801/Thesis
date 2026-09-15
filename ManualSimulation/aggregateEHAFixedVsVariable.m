function summary = aggregateEHAFixedVsVariable()
% aggregateEHAFixedVsVariable.m
% Scans results/EHA_fixedVsVariable/task_*.mat (each written the moment
% its task finishes) and rebuilds summary.csv/.mat -- safe to run at any
% time, including mid-sweep or after an interrupted session. Also prints
% each variant's probability-weighted (expected, across the Humboldt wave
% climate) mechRGP/elecRGP.
%
% Calls: none
% Called by: Run_EHA_FixedVsVariable.m (automatically); also fine to call
%   directly/standalone at any time to check progress.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','EHA_fixedVsVariable');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);
if isempty(files)
    summary = table();
    warning('aggregateEHAFixedVsVariable:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.label}, o.seaStateIdx, o.Hs, o.Tp, o.probability, ...
        o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, o.aveLoss/1e3, ...
        'VariableNames', {'Variant','SeaStateIdx','Hs_m','Tp_s','Probability', ...
        'MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW','AveLoss_kW'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'Variant','SeaStateIdx'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

expectedLabels = {'EHA_fixed_mech','EHA_fixed_elec','EHA_var_mech','EHA_var_elec'};
nSeaStates = max(summary.SeaStateIdx);
fprintf('\n--- Completion check (assuming %d sea states) ---\n', nSeaStates);
for L = expectedLabels
    have = sum(strcmp(summary.Variant, L{1}));
    fprintf('%-16s %d/%d\n', L{1}, have, nSeaStates);
end

fprintf('\n--- Probability-weighted (expected, across the wave climate) ---\n');
for L = expectedLabels
    mask = strcmp(summary.Variant, L{1});
    if ~any(mask), continue; end
    w = summary.Probability(mask);
    w = w / sum(w); % normalize in case only a subset of sea states has landed so far
    expMechRGP = sum(w .* summary.MechRGP(mask));
    expElecRGP = sum(w .* summary.ElecRGP(mask));
    fprintf('%-16s expected mechRGP=%.4f expected elecRGP=%.4f\n', L{1}, expMechRGP, expElecRGP);
end
end
