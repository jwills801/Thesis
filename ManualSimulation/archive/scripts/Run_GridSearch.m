% Run_GridSearch.m
% Big unattended sweep: for each of the 8 Humboldt sea states, find the
% best highPressure for PassivePump and DHD (2/3/4 rails) via a grid
% search (optimization/optimizePressure.m), and just run EHA directly
% (mech- and elec-optimized) since it has no pressure to optimize. Uses
% each drivetrain family's cylinder area from results/sizedAreas.mat
% (produced by optimization/sizeCylinderArea.m -- run that first).
%
% Designed for an unattended, possibly-interrupted cluster run:
%   - Every (drivetrain, sea state) task is independent -> flattened into
%     one task list and run with parfor (no nested parfor; each task's
%     own pressure grid search stays a sequential for-loop inside
%     optimization/optimizePressure.m).
%   - Each task saves its OWN result file (results/gridSearch/task_NN.mat)
%     the moment it finishes, from inside the parfor body -- NOT only at
%     the end. parfor itself only returns results when the entire loop
%     completes, so if the job is killed partway through, relying on that
%     return value would lose everything computed so far. Per-task files
%     survive regardless of when the job stops.
%   - aggregateGridSearchResults.m (separate, rerunnable file) scans
%     whatever task_*.mat files exist and rebuilds the summary
%     table/CSV -- safe to run at any time, including mid-sweep, to check
%     partial progress.
%
% Calls: wave/humboldtSeaStates.m, optimization/buildDenseSwitchMap.m,
%   optimization/optimizePressure.m, parameters/getParameters.m,
%   wave/generateExcitingTorque.m, control/getControl.m,
%   dynamics/timeLoop.m, evaluation/evaluate.m, aggregateGridSearchResults.m
% Called by: none (top-level entry point for the cluster run)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','gridSearch');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

%% Sea states
seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

%% Sized cylinder areas (must already exist -- run sizeCylinderArea.m first)
sizedAreasFile = fullfile(repoRoot,'results','sizedAreas.mat');
if ~isfile(sizedAreasFile)
    error('Run_GridSearch:missingSizedAreas', ...
        ['results/sizedAreas.mat not found. Run the cylinder-area sizing ' ...
         'step first (see optimization/sizeCylinderArea.m) and save its ' ...
         'output struct (fields: PassivePump, EHA, DHD2, DHD3, DHD4, each ' ...
         'with .capArea/.rodArea) to that file before running this script.']);
end
S = load(sizedAreasFile);
sizedAreas = S.sizedAreas;

%% Pressure grid (Pa) -- kept coarse per your call; bump density here if
% you have spare core-hours later.
pressureGrid = linspace(5e6, 35e6, 6);

%% One dense switch-loss map per DHD rail count, built ONCE and reused
% across every sea state and pressure candidate for that family (see
% optimization/buildDenseSwitchMap.m -- avoids regenerating the expensive
% switch-loss map at every grid point). Cached to disk and reloaded if
% present -- each one takes ~28 min to build, so a retry after a crash
% (see the auto-retry loop in run_overnight.sh) must not rebuild these.
denseSwitchMapsFile = fullfile(repoRoot,'results','denseSwitchMaps.mat');
if isfile(denseSwitchMapsFile)
    fprintf('Loading cached dense switch maps from %s\n', denseSwitchMapsFile);
    D = load(denseSwitchMapsFile);
    denseSwitchMaps = D.denseSwitchMaps;
else
denseSwitchMaps = struct();
for nRails = [2 3 4]
    fname = sprintf('DHD%d', nRails);
    hydTmp = struct('capArea', sizedAreas.(fname).capArea, ...
                     'rodArea', sizedAreas.(fname).rodArea, 'stroke', 5);
    fprintf('Building dense switch map for %s (capArea=%.5f)...\n', fname, hydTmp.capArea);
    denseSwitchMaps.(fname) = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10);
end
save(denseSwitchMapsFile,'denseSwitchMaps'); % reused by the A* horizon study (separate process) and by a retry after a crash
end

%% Build the flat task list (one row per drivetrain-family x sea-state)
tasks = struct('label',{}, 'drive',{}, 'controller',{}, 'considerLosses',{}, ...
    'pressure_rails',{}, 'capArea',{}, 'rodArea',{}, 'ehaFixedDisplacement',{}, ...
    'shaftInertia',{}, 'needsPressureOpt',{}, 'denseSwitchMap',{}, ...
    'seaStateIdx',{}, 'Hs',{}, 'Tp',{}, 'probability',{});

familyDefs = {
    'PassivePump', 'PassivePump', 'CoulombDamping', NaN, 2,   true,  []
    'DHD2',        'DHD',         'MPC_Astar',      1,   2,   true,  denseSwitchMaps.DHD2
    'DHD3',        'DHD',         'MPC_Astar',      1,   3,   true,  denseSwitchMaps.DHD3
    'DHD4',        'DHD',         'MPC_Astar',      1,   4,   true,  denseSwitchMaps.DHD4
    'EHA_mech',    'EHA',         'MPC_QP',         0,   NaN, false, []
    'EHA_elec',    'EHA',         'MPC_QP',         1,   NaN, false, []
    };

for fIdx = 1:size(familyDefs,1)
    [label, drive, controller, considerLosses, pressure_rails, needsPressureOpt, denseMap] = familyDefs{fIdx,:};
    switch label
        case 'PassivePump'
            areaKey = 'PassivePump';
        case {'DHD2','DHD3','DHD4'}
            areaKey = label;
        case {'EHA_mech','EHA_elec'}
            areaKey = 'EHA';
    end
    capArea = sizedAreas.(areaKey).capArea;
    rodArea = sizedAreas.(areaKey).rodArea;

    for s = 1:nSeaStates
        t = struct();
        t.label = label; t.drive = drive; t.controller = controller;
        t.considerLosses = considerLosses; t.pressure_rails = pressure_rails;
        t.capArea = capArea; t.rodArea = rodArea;
        t.ehaFixedDisplacement = strcmp(drive,'EHA');
        t.shaftInertia = 5; % kg*m^2, tex's assumption -- see diagnostics/ReadMe.md
        t.needsPressureOpt = needsPressureOpt;
        t.denseSwitchMap = denseMap;
        t.seaStateIdx = s; t.Hs = seaStates.Hs(s); t.Tp = seaStates.Tp(s);
        t.probability = seaStates.probability(s);
        tasks(end+1) = t; %#ok<AGROW>
    end
end

nTasks = numel(tasks);
fprintf('Built %d tasks (%d families x %d sea states). Starting parfor.\n', nTasks, size(familyDefs,1), nSeaStates);

%% Ensure a parallel pool exists, IF Parallel Computing Toolbox is
% available. parfor below runs fine either way -- with no toolbox (this
% machine currently has none installed) it silently falls back to plain
% sequential execution, which is exactly what we want here anyway: this
% machine's 2 concurrent MATLAB processes have deadlocked before
% (license/service-host contention, see earlier in this session), so
% running strictly single-process is the safe choice locally regardless
% of core count. On MSI (interactive-long, 2 cores, toolbox available),
% set nWorkers=2 to get real parallelism -- that's the only line that
% needs to change between environments.
nWorkers = 4; % only used if Parallel Computing Toolbox is actually available
% NOTE: gcp() itself throws a hard error when the toolbox isn't
% installed (license('test',...) checks license entitlement, NOT
% whether the toolbox is actually installed -- it can return true while
% gcp/parpool still error). So the whole thing needs to be inside one
% try/catch, not gated by the license check.
try
    if isempty(gcp('nocreate'))
        parpool('local', nWorkers);
    end
catch ME
    fprintf('No parallel pool available (%s) -- running sequentially (parfor falls back automatically).\n', ME.message);
end

%% Run every task, saving each one's result immediately on completion
parfor i = 1:nTasks
    task = tasks(i); %#ok<PFBNS>
    outFile = fullfile(resultsDir, sprintf('task_%03d.mat', i));
    if isfile(outFile)
        continue % already done (e.g. resuming after an earlier interrupted run)
    end

    runParams = struct('drive',task.drive,'controller',task.controller, ...
        'capArea',task.capArea,'rodArea',task.rodArea);
    if ~isnan(task.considerLosses)
        runParams.considerLosses = task.considerLosses;
    end
    if ~isnan(task.pressure_rails)
        runParams.pressure_rails = task.pressure_rails;
    end
    if task.ehaFixedDisplacement
        runParams.ehaFixedDisplacement = true;
        runParams.shaftInertia = task.shaftInertia;
    end
    if strcmp(task.drive,'DHD')
        % Cap the A* search's per-control-window iteration budget --
        % see control/MPC_Astar.m's comment. Without this, a single
        % DHD4 (nU=16) simulation was observed taking many hours (the
        % node-list re-sort is O(n^2) and n grows with iterMax*(nU-1)).
        % This bounds worst-case runtime; DHD results from this run
        % reflect this tighter search budget, not the algorithm's full
        % potential -- flagged in the final summary too.
        runParams.astarIterMax = 20; % tested: ~5.4min/full-sim for DHD4 vs "many hours, possibly days" unbounded
    end

    seaState = struct('Hs',task.Hs,'Tp',task.Tp);
    if isfield(runParams,'astarIterMax')
        astarIterMaxUsed = runParams.astarIterMax;
    else
        astarIterMaxUsed = NaN; % not applicable (PassivePump/EHA don't use MPC_Astar)
    end
    out = struct('label',task.label,'seaStateIdx',task.seaStateIdx, ...
        'Hs',task.Hs,'Tp',task.Tp,'probability',task.probability,'astarIterMax',astarIterMaxUsed);

    if task.needsPressureOpt
        runParams.highPressure = 35e6; % overwritten per grid point inside optimizePressure
        res = optimizePressure(runParams, seaState, pressureGrid, struct(), task.denseSwitchMap);
        out.bestPressure = res.bestPressure;
        out.mechRGP = res.bestMechRGP;
        out.elecRGP = res.bestElecRGP;
        out.aveMechPow = res.bestAveMechPow;
        out.aveElecPow = res.bestAveElecPow;
        out.nAstarCapHits = res.totalAstarCapHits; % 0 for PassivePump (no A* search)
        out.fullGrid = res; % keep the whole curve too, for later inspection
    else
        params = getParameters(runParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = task.Hs;
        params.simu.peakPeriod = task.Tp;
        wave = generateExcitingTorque(params);
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);
        ev = evaluate(params,dyn,ctrl);
        out.bestPressure = NaN; % not applicable
        out.mechRGP = ev.mechRGP;
        out.elecRGP = ev.elecRGP;
        out.aveMechPow = ev.aveMechPow;
        out.aveElecPow = ev.aveElecPow;
        out.nAstarCapHits = ev.nAstarCapHits; % 0 -- EHA doesn't use MPC_Astar
    end

    parsave(outFile, out);
    fprintf('[%d/%d] %s, sea state %d (Hs=%.2f,Tp=%.2f): mechRGP=%.3f elecRGP=%.3f\n', ...
        i, nTasks, task.label, task.seaStateIdx, task.Hs, task.Tp, out.mechRGP, out.elecRGP);
end

fprintf('All tasks submitted/completed. Aggregating...\n');
aggregateGridSearchResults();
plotPressureSweeps();

fprintf('\nMain sweep done. The A* horizon-length study runs as a separate\n');
fprintf('step (run_overnight.sh) -- each horizon length is its own OS\n');
fprintf('process so a slow one can be killed with a real timeout, which\n');
fprintf('this machine''s lack of Parallel Computing Toolbox rules out doing\n');
fprintf('from inside a single MATLAB session (see prepareAstarTasks.m /\n');
fprintf('runSingleAstarTask.m and run_overnight.sh''s header comments).\n');

function parsave(outFile, out) %#ok<INUSD>
% parfor bodies can't use "save" directly (it looks at the caller's
% workspace, which doesn't work across parallel workers) -- this thin
% wrapper does the same immediate-per-task-save job explicitly instead.
save(outFile,'out');
end
