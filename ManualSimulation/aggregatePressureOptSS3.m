function summary = aggregatePressureOptSS3()
% aggregatePressureOptSS3.m
% Scans results/pressureOptSS3/task_*.mat (each written the moment its
% task finishes) and rebuilds summary.csv/.mat -- safe to run at any
% time, including mid-sweep or after an interrupted session.
%
% Calls: none
% Called by: Run_PressureOptSS3.m (automatically); also fine to call
%   directly/standalone at any time to check progress.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','pressureOptSS3');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);
if isempty(files)
    summary = table();
    warning('aggregatePressureOptSS3:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.family}, o.rails, o.pressure/1e6, o.mAstar, o.mechRGP, o.elecRGP, ...
        o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.elapsedSec, ...
        'VariableNames', {'Family','Rails','Pressure_MPa','MAstar','MechRGP','ElecRGP', ...
        'AveMechPow_kW','AveElecPow_kW','AstarCapHits','ElapsedSec'});
end
summary = sortrows(vertcat(rows{:}), {'Family','Pressure_MPa'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

families = unique(summary.Family);
fprintf('\n--- Completion check (expect 11 pressure points each) ---\n');
for f = families'
    have = sum(strcmp(summary.Family, f{1}));
    fprintf('%-6s %d/11\n', f{1}, have);
end

fprintf('\n--- Best pressure per family (by elecRGP) ---\n');
for f = families'
    rows2 = strcmp(summary.Family, f{1});
    Tf = summary(rows2,:);
    [bestRGP, idx] = max(Tf.ElecRGP);
    fprintf('%-6s best: %.1fMPa -> elecRGP=%.4f mechRGP=%.4f  [%d/11 points completed]\n', ...
        f{1}, Tf.Pressure_MPa(idx), bestRGP, Tf.MechRGP(idx), height(Tf));
end

if any(summary.AstarCapHits > 0)
    fprintf('\nWARNING: some tasks hit the A* iteration cap -- check which and investigate.\n');
end
end
