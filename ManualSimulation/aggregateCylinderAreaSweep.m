function summary = aggregateCylinderAreaSweep()
% aggregateCylinderAreaSweep.m
% Scans results/cylinderAreaSweep/task_*.mat (each written the moment
% its task finishes) and rebuilds summary.csv/.mat -- safe to run at any
% time, including mid-sweep or after an interrupted session.
%
% Calls: none
% Called by: Run_CylinderAreaSweep.m (automatically); also fine to call
%   directly/standalone at any time.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','cylinderAreaSweep');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);
if isempty(files)
    summary = table();
    warning('aggregateCylinderAreaSweep:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.family}, o.rails, o.diamIn, o.capArea, o.mAstar, o.mechRGP, o.elecRGP, ...
        o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.elapsedSec, ...
        'VariableNames', {'Family','Rails','DiamIn','CapArea_m2','MAstar','MechRGP','ElecRGP', ...
        'AveMechPow_kW','AveElecPow_kW','AstarCapHits','ElapsedSec'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'Family','DiamIn'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

families = unique(summary.Family);
fprintf('\n--- Best bore diameter per family (by elecRGP) ---\n');
for f = families'
    rows = strcmp(summary.Family, f{1});
    Tf = summary(rows,:);
    [bestRGP, idx] = max(Tf.ElecRGP);
    fprintf('%-6s best: %.1fin (capArea=%.6f) -> elecRGP=%.4f  [%d/%d diameters completed]\n', ...
        f{1}, Tf.DiamIn(idx), Tf.CapArea_m2(idx), bestRGP, height(Tf), 10);
end
end
