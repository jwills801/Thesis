function [astarSummary, pressureSummary] = aggregateMAstarConvAndPressureOpt()
% aggregateMAstarConvAndPressureOpt.m
% Scans results/mAstarConvAndPressureOpt/task_*.mat (each written the
% moment its task finishes) and rebuilds both summary tables -- safe to
% run at any time, including mid-sweep or after an interrupted session.
%
% Calls: none
% Called by: Run_MAstarConvAndPressureOpt.m (automatically); also fine to
%   call directly/standalone at any time to check progress.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','mAstarConvAndPressureOpt');

%% Part A: m_Astar convergence
files = dir(fullfile(resultsDir,'task_astarConv_*.mat'));
fprintf('Found %d astarConv task file(s)\n', numel(files));
if isempty(files)
    astarSummary = table();
else
    rows = cell(numel(files),1);
    for i = 1:numel(files)
        d = load(fullfile(files(i).folder,files(i).name));
        o = d.out;
        rows{i} = table({o.family}, o.rails, o.mAstar, o.mechRGP, o.elecRGP, ...
            o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.elapsedSec, ...
            'VariableNames', {'Family','Rails','MAstar','MechRGP','ElecRGP', ...
            'AveMechPow_kW','AveElecPow_kW','AstarCapHits','ElapsedSec'});
    end
    astarSummary = sortrows(vertcat(rows{:}), {'Family','MAstar'});
    writetable(astarSummary, fullfile(resultsDir,'astarConvSummary.csv'));
    save(fullfile(resultsDir,'astarConvSummary.mat'),'astarSummary');
    fprintf('Wrote astarConvSummary.csv/.mat\n');
    disp(astarSummary);

    fprintf('\n--- astarConv completion check (expect 10 mAstar values each) ---\n');
    for fam = unique(astarSummary.Family)'
        rows2 = strcmp(astarSummary.Family, fam{1});
        fprintf('%-6s mAstar values done: %s\n', fam{1}, mat2str(sort(astarSummary.MAstar(rows2))'));
    end
end

%% Part B: pressure optimization
files = dir(fullfile(resultsDir,'task_pressureOpt_*.mat'));
fprintf('\nFound %d pressureOpt task file(s)\n', numel(files));
if isempty(files)
    pressureSummary = table();
    return
end
rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table({o.family}, o.rails, o.seaStateIdx, o.Hs, o.Tp, o.probability, ...
        o.bestPressure/1e6, o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, ...
        'VariableNames', {'Family','Rails','SeaStateIdx','Hs_m','Tp_s','Probability', ...
        'BestPressure_MPa','MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW','AstarCapHits'});
end
pressureSummary = sortrows(vertcat(rows{:}), {'Family','SeaStateIdx'});
writetable(pressureSummary, fullfile(resultsDir,'pressureOptSummary.csv'));
save(fullfile(resultsDir,'pressureOptSummary.mat'),'pressureSummary');
fprintf('Wrote pressureOptSummary.csv/.mat\n');
disp(pressureSummary);

fprintf('\n--- pressureOpt completion check (expect 8 sea states each) ---\n');
families = unique(pressureSummary.Family);
for fam = families'
    have = sum(strcmp(pressureSummary.Family, fam{1}));
    fprintf('%-12s %d/8\n', fam{1}, have);
end

fprintf('\n--- Probability-weighted (expected, across the wave climate) ---\n');
for fam = families'
    mask = strcmp(pressureSummary.Family, fam{1});
    if ~any(mask), continue; end
    w = pressureSummary.Probability(mask);
    w = w / sum(w);
    expMechRGP = sum(w .* pressureSummary.MechRGP(mask));
    expElecRGP = sum(w .* pressureSummary.ElecRGP(mask));
    fprintf('%-12s expected mechRGP=%.4f expected elecRGP=%.4f\n', fam{1}, expMechRGP, expElecRGP);
end

if any(pressureSummary.AstarCapHits > 0)
    fprintf('\nWARNING: some DHD pressureOpt tasks hit the A* iteration cap -- check which and investigate.\n');
end
end
