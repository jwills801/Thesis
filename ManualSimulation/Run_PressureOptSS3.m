% Run_PressureOptSS3.m
% Quick pressure optimization at sea state 3 for DHD2/DHD3/DHD4 (each at
% its own newly-locked bore -- PassivePump=13in/DHD2=16in/DHD3=17in/
% DHD4=18in, see updateFamilyAreas.m, results/sizedAreas.mat), mAstar
% fixed at 5. PassivePump deliberately excluded here -- planned for the
% next grid search instead.
%
% Run this BEFORE Run_MAstarConvAndPressureOpt.m: that script's Part A
% (mAstar convergence at SS3) assumes a fixed 35MPa, which was chosen for
% SS8 (a much more energetic sea state) -- likely not right for SS3.
% This finds SS3's own near-optimal pressure per family so the
% convergence study can use a sensible value instead.
%
% Flattened to one task per (family, pressure) point -- not one task per
% family running an 11-point grid internally -- specifically to make use
% of a full 32-core allocation: 3 families x 11 pressure points = 33
% tasks (a finer grid than the usual 6-point convention, since this is a
% cheap single-sea-state check with cores to spare).
%
% Builds (or reuses, if already built) the same shared switch-loss map
% used by Run_MAstarConvAndPressureOpt.m (results/denseSwitchMap_shared_
% 100ms.mat, sized to DHD4's 18in bore) -- running this script first
% means that map is already cached when the later script runs.
%
% Writes results/pressureOptSS3/task_<family>_p<pressureMPa>MPa.mat (one
% per task, saved immediately) and, via aggregatePressureOptSS3.m,
% results/pressureOptSS3/summary.csv.
%
% Calls: optimization/buildDenseSwitchMap.m, aggregatePressureOptSS3.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','pressureOptSS3');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

seaStates = humboldtSeaStates();
ss3Idx = 3;
Hs3 = seaStates.Hs(ss3Idx); Tp3 = seaStates.Tp(ss3Idx);

CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;
MASTAR_FIXED = 5;
pressureGrid = linspace(5e6, 35e6, 11); % finer than the usual 6-point grid -- cheap here, cores to spare

dhdFamilies = {'DHD2',2; 'DHD3',3; 'DHD4',4};

%% Build (or load) ONE shared switch-loss map, sized to the largest DHD bore (18in, DHD4)
mapFile = fullfile(repoRoot,'results','denseSwitchMap_shared_100ms.mat');
if isfile(mapFile)
    fprintf('Loading cached shared switch map from %s\n', mapFile);
    D = load(mapFile);
    sharedMap = D.denseMap;
else
    fprintf('Building 100ms dense switch map sized to DHD4''s 18in bore (capArea=%.6f)...\n', sizedAreas.DHD4.capArea);
    tic;
    hydTmp = struct('capArea',sizedAreas.DHD4.capArea,'rodArea',sizedAreas.DHD4.rodArea,'stroke',5);
    denseMap = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.1); %#ok<NASGU>
    fprintf('Build took %.1f minutes\n', toc/60);
    save(mapFile,'denseMap');
    sharedMap = denseMap;
end

%% Flatten into one task per (family, pressure)
tasks = struct('family',{},'rails',{},'pressure',{});
for f = 1:size(dhdFamilies,1)
    for p = pressureGrid
        tasks(end+1) = struct('family',dhdFamilies{f,1},'rails',dhdFamilies{f,2},'pressure',p); %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d tasks (%d families x %d pressure points) at sea state 3 (Hs=%.2f, Tp=%.2f).\n', ...
    nTasks, size(dhdFamilies,1), numel(pressureGrid), Hs3, Tp3);

nWorkers = min(nTasks, 32);
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
    outFile = fullfile(resultsDir, sprintf('task_%s_p%.1fMPa.mat', task.family, task.pressure/1e6));
    if isfile(outFile), continue; end

    runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
        'considerLosses',1,'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea, ...
        'highPressure',task.pressure,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',MASTAR_FIXED);
    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs3;
    params.simu.peakPeriod = Tp3;
    params.hyd.switchMap = sharedMap; %#ok<PFBNS>

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    tic;
    dyn = timeLoop(params,wave,ctrl);
    elapsed = toc;
    ev = evaluate(params,dyn,ctrl);

    out = struct('family',task.family,'rails',task.rails,'pressure',task.pressure,'mAstar',MASTAR_FIXED, ...
        'Hs',Hs3,'Tp',Tp3,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP,'aveMechPow',ev.aveMechPow, ...
        'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits,'elapsedSec',elapsed);
    parsave(outFile, out);
    fprintf('[%d/%d] %s p=%.1fMPa: mechRGP=%.3f elecRGP=%.3f capHits=%d (%.1f sec)\n', ...
        i, nTasks, task.family, task.pressure/1e6, out.mechRGP, out.elecRGP, out.nAstarCapHits, elapsed);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregatePressureOptSS3();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
