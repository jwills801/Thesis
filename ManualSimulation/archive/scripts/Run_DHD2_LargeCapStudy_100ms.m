% Run_DHD2_LargeCapStudy_100ms.m
% Same as Run_DHD2_LargeCapStudy.m (DHD2 pressure grid search, A*
% iteration cap raised to 100000) but with the coarse control step
% dropped from 200ms to 100ms (runParams.controlDT), to test whether
% DHD2 trailing PassivePump is partly a control-cadence effect -- a
% rail decision can currently only change every 200ms, versus
% PassivePump's effectively-continuous response.
%
% Builds (or reuses, if already cached) its OWN switch-loss map at
% switchTime=0.1s -- the 200ms dense map (results/denseSwitchMaps.mat)
% is NOT valid here, see models/makeSwitchLossMap.m's comment on why
% finalTime must match ctrl.timeHorizon. Cached separately to
% results/denseSwitchMaps_100ms.mat so both the 100ms and 200ms dense
% maps stay available without rebuilding either.
%
% Writes to results/gridSearchLargeCap_100ms/ (kept separate from the
% 200ms study's results/gridSearchLargeCap/ for a direct comparison).
%
% Calls: wave/humboldtSeaStates.m, optimization/buildDenseSwitchMap.m,
%   optimization/optimizePressure.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','gridSearchLargeCap_100ms');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

CONTROL_DT = 0.1;
M_ASTAR = 10; % 10*0.1s = 1.0s of foresight, matching the 200ms study's 5*0.2s -- holds
              % lookahead TIME constant so this isolates switching cadence, not foresight window

denseMapFile = fullfile(repoRoot,'results','denseSwitchMaps_100ms.mat');
if isfile(denseMapFile)
    fprintf('Loading cached 100ms dense switch map from %s\n', denseMapFile);
    D = load(denseMapFile);
    denseMapDHD2 = D.denseMapDHD2;
else
    fprintf('Building 100ms dense switch map for DHD2 (capArea=%.5f)...\n', sizedAreas.DHD2.capArea);
    hydTmp = struct('capArea', sizedAreas.DHD2.capArea, 'rodArea', sizedAreas.DHD2.rodArea, 'stroke', 5);
    denseMapDHD2 = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, CONTROL_DT);
    save(denseMapFile,'denseMapDHD2');
end

pressureGrid = linspace(5e6, 35e6, 6); % same grid as the 200ms study, for a fair comparison

ASTAR_ITER_MAX = 100000;

nWorkers = min(nSeaStates, 4);
slurmCpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(slurmCpus)
    slurmCpus = str2double(getenv('SLURM_CPUS_ON_NODE'));
end
if ~isnan(slurmCpus) && slurmCpus > 0
    nWorkers = min(nSeaStates, slurmCpus);
end
try
    if isempty(gcp('nocreate'))
        parpool('local', nWorkers);
    end
catch ME
    fprintf('No parallel pool available (%s) -- running sequentially.\n', ME.message);
end

parfor s = 1:nSeaStates
    outFile = fullfile(resultsDir, sprintf('DHD2_seaState_%02d.mat', s));
    if isfile(outFile)
        continue
    end
    runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',2, ...
        'considerLosses',1,'capArea',sizedAreas.DHD2.capArea,'rodArea',sizedAreas.DHD2.rodArea, ...
        'highPressure',35e6,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',M_ASTAR);
    seaState = struct('Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s));
    res = optimizePressure(runParams, seaState, pressureGrid, struct(), denseMapDHD2);
    out = struct('seaStateIdx',s,'Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s), ...
        'probability',seaStates.probability(s),'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT, ...
        'bestPressure',res.bestPressure,'mechRGP',res.bestMechRGP,'elecRGP',res.bestElecRGP, ...
        'aveMechPow',res.bestAveMechPow,'aveElecPow',res.bestAveElecPow, ...
        'nAstarCapHits',res.totalAstarCapHits,'fullGrid',res);
    parsave(outFile, out);
    fprintf('[%d/%d] DHD2 sea state %d (Hs=%.2f,Tp=%.2f): mechRGP=%.3f elecRGP=%.3f capHits=%d\n', ...
        s, nSeaStates, s, seaStates.Hs(s), seaStates.Tp(s), out.mechRGP, out.elecRGP, out.nAstarCapHits);
end

%% Aggregate
files = dir(fullfile(resultsDir,'DHD2_seaState_*.mat'));
rows = cell(numel(files),1);
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    rows{i} = table(o.seaStateIdx, o.Hs, o.Tp, o.probability, o.bestPressure/1e6, ...
        o.mechRGP, o.elecRGP, o.aveMechPow/1e3, o.aveElecPow/1e3, o.nAstarCapHits, o.astarIterMax, ...
        'VariableNames', {'SeaStateIdx','Hs_m','Tp_s','Probability','BestPressure_MPa', ...
        'MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW','AstarCapHits','AstarIterMax'});
end
summary = vertcat(rows{:});
summary = sortrows(summary,'SeaStateIdx');
writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
disp(summary);

totalCapHits = sum(summary.AstarCapHits);
if totalCapHits > 0
    fprintf('\nWARNING: %d total A* cap hits across the DHD2 100ms sweep even at astarIterMax=%d.\n', totalCapHits, ASTAR_ITER_MAX);
else
    fprintf('\nConfirmed: astarIterMax=%d was never hit anywhere in the DHD2 100ms sweep.\n', ASTAR_ITER_MAX);
end

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
