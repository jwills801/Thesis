% Run_CylinderAreaSweep.m
% Cylinder-bore-diameter sweep at sea state 8 (Hs=4.995m, Tp=12.725s),
% run SEPARATELY for PassivePump, DHD2, DHD3, and DHD4 -- i.e. each
% family gets its own independent bore-size optimum, not assumed to
% share one (the earlier convention that all three DHD families use
% identical capArea/rodArea was a sizing-methodology shortcut, not a
% physical requirement; same logic extends to PassivePump). m_Astar
% (DHD only) is fixed at 5, pressure fixed at 35MPa for every family
% (not swept) -- this sweep is purely about bore size. Pressure gets
% re-optimized later, per sea state, at whatever fixed area each family
% lands on here -- that later step is where PassivePump's traditionally
% much lower ~20.6MPa operating point would actually get reconsidered,
% not here.
%
% boreDiametersIn (edit this to change the sweep): 10 points, 10-25in --
% sized so 4 families x 10 diameters = 40 tasks fits comfortably under a
% 32-core allocation without extending wall-clock time much versus the
% original 3-family/30-task version (PassivePump tasks are much cheaper
% than DHD's -- no A* search, no switch-loss map). Brackets both
% previously-tested reference points (14.70in = current DHD sizing,
% 18in = the earlier one-off DHD3/DHD4 test) plus smaller/larger
% candidates on either side.
%
% Only ONE switch-loss map is built for the WHOLE sweep (not one per
% bore diameter, and not one per DHD family) -- makeSwitchLossMap.m's
% actual loss physics depends only on the real flow/volume values at
% each grid point, never on capArea directly; capArea only sets how far
% the velA/vol axes need to extend. Building that one map sized to the
% LARGEST bore diameter in the sweep covers every smaller bore's actual
% (narrower) flow/volume range too, so it's reused across all 10
% diameters and all three DHD families (rail count is a control-layer
% choice, not a cylinder-sizing one -- see
% optimization/buildDenseSwitchMap.m's header). PassivePump needs no
% such map at all: its check-valve loss (hyd.switchMap.valveConstant) is
% sized directly from capArea/vMax inside getHydraulic.m.
%
% Writes results/cylinderAreaSweep/task_<family>_<diamIn>in.mat (one per
% task, saved immediately) and, via aggregateCylinderAreaSweep.m,
% results/cylinderAreaSweep/summary.csv.
%
% Calls: optimization/buildDenseSwitchMap.m, aggregateCylinderAreaSweep.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','cylinderAreaSweep');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

Hs8 = 4.99525316455696;
Tp8 = 12.7246835443038;
PRESSURE = 35e6;
CONTROL_DT = 0.1;
MASTAR = 5;
ASTAR_ITER_MAX = 100000;

boreDiametersIn = [10 11.7 13.3 14.7 16.7 18 20 21.7 23.3 25]; % EDIT HERE to change the sweep
families = {'PassivePump',NaN; 'DHD2',2; 'DHD3',3; 'DHD4',4}; % rails: NaN = not applicable (PassivePump)

in2m = 0.0254;
diamToAreas = @(d_in) deal(pi*(d_in*in2m/2)^2, pi*(d_in*in2m/2)^2/1.5); % [capArea, rodArea], 1.5:1 ratio

%% Build (or load) ONE switch-loss map, sized to the LARGEST bore diameter
% makeSwitchLossMap.m's actual loss physics (valve orifice equations,
% pressure buildup, valveConstant) depends only on the real flow (velA,
% m^3/s) and real volume (vol, m^3) at each grid point, plus the fixed
% pressure rails and fluid properties -- capArea never enters the
% physics directly, only the *range* of the velA/vol axes (+-1.5*capArea
% and hoseVolume+stroke*capArea respectively). So one map built wide
% enough to cover the largest bore's flow/volume range is valid for
% every smaller bore too (interpn just uses a subset of that range) --
% no need to rebuild per diameter.
maxDiamIn = max(boreDiametersIn);
[capAreaMax, rodAreaMax] = diamToAreas(maxDiamIn);
mapFile = fullfile(repoRoot,'results',sprintf('denseSwitchMap_%gin_100ms.mat', maxDiamIn));
if isfile(mapFile)
    fprintf('Loading cached switch map (sized to largest bore, %gin) from %s\n', maxDiamIn, mapFile);
    D = load(mapFile);
    sharedMap = D.denseMap;
else
    fprintf('Building 100ms dense switch map sized to largest bore (%gin, capArea=%.6f)...\n', maxDiamIn, capAreaMax);
    tic;
    hydTmp = struct('capArea',capAreaMax,'rodArea',rodAreaMax,'stroke',5);
    denseMap = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.1); %#ok<NASGU>
    fprintf('Build took %.1f minutes\n', toc/60);
    save(mapFile,'denseMap');
    sharedMap = denseMap;
end

%% Flatten into one task list (family x bore diameter)
tasks = struct('family',{},'rails',{},'diamIn',{},'capArea',{},'rodArea',{});
for f = 1:size(families,1)
    for i = 1:numel(boreDiametersIn)
        d = boreDiametersIn(i);
        [capArea, rodArea] = diamToAreas(d);
        tasks(end+1) = struct('family',families{f,1},'rails',families{f,2}, ...
            'diamIn',d,'capArea',capArea,'rodArea',rodArea); %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d tasks (%d families x %d bore diameters).\n', nTasks, size(families,1), numel(boreDiametersIn));

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
    outFile = fullfile(resultsDir, sprintf('task_%s_%gin.mat', task.family, task.diamIn));
    if isfile(outFile), continue; end

    if strcmp(task.family,'PassivePump')
        runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
            'capArea',task.capArea,'rodArea',task.rodArea, ...
            'highPressure',PRESSURE,'controlDT',CONTROL_DT);
    else
        runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
            'considerLosses',1,'capArea',task.capArea,'rodArea',task.rodArea, ...
            'highPressure',PRESSURE,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',MASTAR);
    end
    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs8;
    params.simu.peakPeriod = Tp8;
    if ~strcmp(task.family,'PassivePump')
        params.hyd.switchMap = sharedMap; %#ok<PFBNS>
    end

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    tic;
    dyn = timeLoop(params,wave,ctrl);
    elapsed = toc;
    ev = evaluate(params,dyn,ctrl);

    mAstarOut = MASTAR; if strcmp(task.family,'PassivePump'), mAstarOut = NaN; end
    out = struct('family',task.family,'rails',task.rails,'diamIn',task.diamIn, ...
        'capArea',task.capArea,'rodArea',task.rodArea,'mAstar',mAstarOut,'Hs',Hs8,'Tp',Tp8, ...
        'pressure',PRESSURE,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP,'aveMechPow',ev.aveMechPow, ...
        'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits,'elapsedSec',elapsed);
    parsave(outFile, out);
    fprintf('[%d/%d] %s %gin: mechRGP=%.3f elecRGP=%.3f capHits=%d (%.1f sec)\n', ...
        i, nTasks, task.family, task.diamIn, out.mechRGP, out.elecRGP, out.nAstarCapHits, elapsed);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateCylinderAreaSweep();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
