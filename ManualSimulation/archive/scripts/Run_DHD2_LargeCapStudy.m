% Run_DHD2_LargeCapStudy.m
% Re-runs DHD2's pressure grid search (same procedure as Run_GridSearch.m)
% with astarIterMax raised from 20 to 100000 -- a budget confirmed by a
% timing probe (60s sim at the worst-case sea state, 0.122 sec/simulated-
% sec, zero cap hits) to run fast and reflect A*'s natural convergence
% rather than a truncated search. DHD2 goes first because it's the
% cheapest DHD configuration (nU=4); DHD3/DHD4 need their own timing
% probe before committing to a similarly large cap, since MPC_Astar.m's
% O(n^2) node-list resort scales badly with nU (this was the whole
% reason astarIterMax=20 was used for the original Run_GridSearch.m
% sweep -- see its comments).
%
% Writes to results/gridSearchLargeCap/ (a NEW directory, not
% results/gridSearch/) so the astarIterMax=20 baseline is preserved for
% a direct before/after comparison.
%
% Calls: wave/humboldtSeaStates.m, optimization/optimizePressure.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','gridSearchLargeCap');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

D = load(fullfile(repoRoot,'results','denseSwitchMaps.mat'));
denseMapDHD2 = D.denseSwitchMaps.DHD2;

pressureGrid = linspace(5e6, 35e6, 6); % same grid as the capped baseline, for a fair comparison

ASTAR_ITER_MAX = 100000; % confirmed fast (0.122 sec/simulated-sec) and never hit at the worst-case sea state

% Worker count: matches nSeaStates (8) when possible, since this loop's
% only parallelism is one worker per sea state -- more workers than that
% buys nothing here. Reads SLURM_CPUS_PER_TASK (set by srun --cpus-per-
% task) or SLURM_CPUS_ON_NODE (set by an interactive salloc allocation
% instead) when present, else falls back to 4 (this machine's local
% default, unchanged from before). If no Parallel Computing Toolbox is
% available (true for this machine), gcp/parpool throws and the catch
% below falls back to plain serial execution -- parfor degrades
% gracefully either way.
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
        'highPressure',35e6,'astarIterMax',ASTAR_ITER_MAX);
    seaState = struct('Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s));
    res = optimizePressure(runParams, seaState, pressureGrid, struct(), denseMapDHD2);
    out = struct('seaStateIdx',s,'Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s), ...
        'probability',seaStates.probability(s),'astarIterMax',ASTAR_ITER_MAX, ...
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
    fprintf('\nWARNING: %d total A* cap hits across the DHD2 sweep even at astarIterMax=%d -- investigate before trusting these as uncapped results.\n', totalCapHits, ASTAR_ITER_MAX);
else
    fprintf('\nConfirmed: astarIterMax=%d was never hit anywhere in the DHD2 sweep -- these results reflect A*''s natural convergence, not a truncated search.\n', ASTAR_ITER_MAX);
end

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
