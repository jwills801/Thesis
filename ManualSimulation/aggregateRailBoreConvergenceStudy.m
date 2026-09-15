function summary = aggregateRailBoreConvergenceStudy()
% aggregateRailBoreConvergenceStudy.m
% Scans results/railBoreConvergence/task_*.mat (each written the moment
% its task finishes) and rebuilds summary.csv/.mat -- safe to run at any
% time, including mid-sweep or after an interrupted session.
%
% Calls: none
% Called by: Run_RailBoreConvergenceStudy.m (automatically); also fine
%   to call directly/standalone at any time.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','railBoreConvergence');

files = dir(fullfile(resultsDir,'task_*.mat'));
fprintf('Found %d completed task file(s) in %s\n', numel(files), resultsDir);
if isempty(files)
    summary = table();
    warning('aggregateRailBoreConvergenceStudy:noTasks','No task_*.mat files found yet.');
    return
end

rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.family}, o.rails, {o.bore}, o.capArea, o.mAstar, o.mechRGP, o.elecRGP, ...
        o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.elapsedSec, ...
        'VariableNames', {'Family','Rails','Bore','CapArea_m2','MAstar','MechRGP','ElecRGP', ...
        'AveMechPow_kW','AveElecPow_kW','AstarCapHits','ElapsedSec'});
end
summary = vertcat(rows{:});
summary = sortrows(summary, {'Family','Bore','MAstar'});

writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('Wrote %s and %s\n', fullfile(resultsDir,'summary.csv'), fullfile(resultsDir,'summary.mat'));
disp(summary);

families = unique(summary.Family);
bores = unique(summary.Bore);
fprintf('\n--- Completion check (expect %d mAstar values each) ---\n', 10);
for f = families'
    for b = bores'
        rows = strcmp(summary.Family,f{1}) & strcmp(summary.Bore,b{1});
        have = sort(summary.MAstar(rows))';
        fprintf('%-6s %-8s mAstar values done: %s\n', f{1}, b{1}, mat2str(have));
    end
end
end
