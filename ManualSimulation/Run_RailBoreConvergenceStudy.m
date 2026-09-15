% Run_RailBoreConvergenceStudy.m
% A* search-horizon (m_Astar) convergence sweep at sea state 8
% (Hs=4.995m, Tp=12.725s -- the sizing/"most energetic" sea state), for
% ALL THREE DHD rail counts (2/3/4) at BOTH the current bore (14.70in,
% from results/sizedAreas.mat) and an 18in bore (capArea=0.164173,
% rodArea=0.109449 m^2). Pressure is fixed at 35MPa throughout (not
% swept -- this is about the horizon-length/bore-size effect, not
% pressure optimization).
%
% Only TWO switch-loss maps are needed (one per bore), not six -- DHD2/
% DHD3/DHD4 share identical capArea/rodArea for a given bore (rail count
% is a control-layer choice, not a cylinder-sizing one; see
% optimization/buildDenseSwitchMap.m's header), so each bore's map is
% built once and reused across all three rail counts. Built fresh here
% if not already cached (~15-20 min each) -- safe to run standalone
% anywhere, including fresh on a cluster with nothing pre-built.
%
% Writes results/railBoreConvergence/task_<family>_<bore>_m<mAstar>.mat
% (one per task, saved immediately) and, via
% aggregateRailBoreConvergenceStudy.m, results/railBoreConvergence/summary.csv.
%
% Calls: optimization/buildDenseSwitchMap.m, aggregateRailBoreConvergenceStudy.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','railBoreConvergence');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

Hs8 = 4.99525316455696;
Tp8 = 12.7246835443038;
PRESSURE = 35e6;
CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;
mAstarValues = [2 4 6 8 10 12 14 16 18 20];

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));

% ONE (capArea,rodArea) pair per bore, shared across DHD2/DHD3/DHD4
bores = struct( ...
    'current', struct('capArea',S.sizedAreas.DHD3.capArea, 'rodArea',S.sizedAreas.DHD3.rodArea, ...
                       'mapFile','denseSwitchMaps_100ms.mat', 'mapVar','denseMapDHD2'), ...
    'in18',    struct('capArea',0.164173, 'rodArea',0.109449, ...
                       'mapFile','denseSwitchMap_18in_100ms.mat', 'mapVar','denseMap18') ...
    );
boreNames = fieldnames(bores);

families = {'DHD2',2; 'DHD3',3; 'DHD4',4};

%% Load or build each bore's 100ms dense switch map (once, before the sweep)
maps = struct();
for b = 1:numel(boreNames)
    name = boreNames{b};
    bore = bores.(name);
    mapPath = fullfile(repoRoot,'results',bore.mapFile);
    if isfile(mapPath)
        fprintf('Loading cached switch map for %s from %s\n', name, mapPath);
        D = load(mapPath);
        maps.(name) = D.(bore.mapVar);
    else
        fprintf('Building 100ms dense switch map for %s (capArea=%.6f)...\n', name, bore.capArea);
        tic;
        hydTmp = struct('capArea',bore.capArea,'rodArea',bore.rodArea,'stroke',5);
        maps.(name) = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.1);
        fprintf('Build took %.1f minutes\n', toc/60);
        s = struct(bore.mapVar, maps.(name));
        save(mapPath, '-struct', 's');
    end
end

%% Flatten into one task list (family x bore x mAstar)
tasks = struct('family',{},'rails',{},'bore',{},'mAstar',{});
for f = 1:size(families,1)
    for b = 1:numel(boreNames)
        for m = mAstarValues
            tasks(end+1) = struct('family',families{f,1},'rails',families{f,2}, ...
                'bore',boreNames{b},'mAstar',m); %#ok<AGROW>
        end
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d tasks (%d families x %d bores x %d mAstar values).\n', ...
    nTasks, size(families,1), numel(boreNames), numel(mAstarValues));

nWorkers = min(nTasks, 4);
slurmCpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(slurmCpus), slurmCpus = str2double(getenv('SLURM_CPUS_ON_NODE')); end
if ~isnan(slurmCpus) && slurmCpus > 0, nWorkers = min(nTasks, slurmCpus); end
try
    if isempty(gcp('nocreate')), parpool('local', nWorkers); end
catch ME
    fprintf('No parallel pool available (%s) -- running sequentially.\n', ME.message);
end

parfor i = 1:nTasks
    task = tasks(i); %#ok<PFBNS>
    outFile = fullfile(resultsDir, sprintf('task_%s_%s_m%02d.mat', task.family, task.bore, task.mAstar));
    if isfile(outFile), continue; end

    bore = bores.(task.bore); %#ok<PFBNS>
    runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
        'considerLosses',1,'capArea',bore.capArea,'rodArea',bore.rodArea, ...
        'highPressure',PRESSURE,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',task.mAstar);
    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs8;
    params.simu.peakPeriod = Tp8;
    params.hyd.switchMap = maps.(task.bore); %#ok<PFBNS>

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    tic;
    dyn = timeLoop(params,wave,ctrl);
    elapsed = toc;
    ev = evaluate(params,dyn,ctrl);

    out = struct('family',task.family,'rails',task.rails,'bore',task.bore,'capArea',bore.capArea, ...
        'mAstar',task.mAstar,'Hs',Hs8,'Tp',Tp8,'pressure',PRESSURE,'mechRGP',ev.mechRGP, ...
        'elecRGP',ev.elecRGP,'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow, ...
        'nAstarCapHits',ev.nAstarCapHits,'elapsedSec',elapsed);
    parsave(outFile, out);
    fprintf('[%d/%d] %s/%s m=%d: mechRGP=%.3f elecRGP=%.3f capHits=%d (%.1f sec)\n', ...
        i, nTasks, task.family, task.bore, task.mAstar, out.mechRGP, out.elecRGP, out.nAstarCapHits, elapsed);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateRailBoreConvergenceStudy();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
