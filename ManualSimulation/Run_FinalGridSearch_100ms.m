% Run_FinalGridSearch_100ms.m
% The "final" pressure grid search: all 8 Humboldt sea states x
% {PassivePump, DHD2, DHD3, DHD4, EHA_mech, EHA_elec}, with the
% getMotorLoss/valve-frequency bug fixes in place, at controlDT=0.1
% (100ms, the chosen coarse time step going forward) and
% astarIterMax=100000 (confirmed to never actually bind at reasonable
% mAstar). Each DHD family uses its OWN mAstar, auto-selected as whichever
% value gave the best elecRGP in results/astarConvergence100ms/summary.csv
% (the sea-state-7 convergence study) -- see that script's header for why
% a single-sea-state study is used for all 8 here.
%
% Requires, all already current as of this script's first run:
%   - results/sizedAreas.mat
%   - results/denseSwitchMaps_100ms.mat (rebuilt with corrected valve
%     dynamics by rebuildSwitchMaps.m -- ONE map, valid for DHD2/3/4 alike)
%   - results/astarConvergence100ms/summary.csv (from
%     Run_AstarConvergenceStudy_100ms.m, run to completion first)
%
% Writes results/gridSearchFinal/task_%03d.mat (one per task, saved the
% moment it finishes) and, via aggregateFinalGridSearchResults.m,
% results/gridSearchFinal/summary.csv.
%
% Calls: wave/humboldtSeaStates.m, optimization/optimizePressure.m,
%   aggregateFinalGridSearchResults.m, plotFinalPressureSweeps.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','gridSearchFinal');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

D = load(fullfile(repoRoot,'results','denseSwitchMaps_100ms.mat'));
denseMap = D.denseMapDHD2; % valid for DHD2/DHD3/DHD4 alike -- see its header comment

%% Each DHD family's mAstar -- MANUALLY chosen after reviewing
% results/astarConvergence100ms/plots/astarConvergence.png (and/or its
% summary.csv) by eye, not auto-picked. A blind argmax on elecRGP risks
% latching onto search noise rather than a genuine plateau -- this
% search has already shown real run-to-run variation on identical inputs
% across hardware (see diagnostics/ReadMe.md). Fill these in once you've
% looked at the convergence plot; leaving any of them empty ([]) errors
% below rather than silently guessing.
MASTAR_CHOICE = struct('DHD2', [], 'DHD3', [], 'DHD4', []);

convFile = fullfile(repoRoot,'results','astarConvergence100ms','summary.csv');
dhdFamilies = {'DHD2','DHD3','DHD4'};
mAstarChoice = struct();
for f = 1:numel(dhdFamilies)
    fam = dhdFamilies{f};
    if isempty(MASTAR_CHOICE.(fam))
        error('Run_FinalGridSearch_100ms:noMAstarChosen', ...
            ['MASTAR_CHOICE.%s is not set -- review %s (and its plot) and fill in a ' ...
             'value at the top of this script before running.'], fam, convFile);
    end
    mAstarChoice.(fam) = MASTAR_CHOICE.(fam);
    fprintf('%s: using manually-chosen mAstar=%d\n', fam, mAstarChoice.(fam));
end

%% Pressure grid (unchanged from the earlier sweep)
pressureGrid = linspace(5e6, 35e6, 6);

CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;

%% Build the flat task list
familyDefs = {
    'PassivePump', 'PassivePump', 'CoulombDamping', NaN, 2,   true,  NaN
    'DHD2',        'DHD',         'MPC_Astar',      1,   2,   true,  mAstarChoice.DHD2
    'DHD3',        'DHD',         'MPC_Astar',      1,   3,   true,  mAstarChoice.DHD3
    'DHD4',        'DHD',         'MPC_Astar',      1,   4,   true,  mAstarChoice.DHD4
    'EHA_mech',    'EHA',         'MPC_QP',         0,   NaN, false, NaN
    'EHA_elec',    'EHA',         'MPC_QP',         1,   NaN, false, NaN
    };

tasks = struct('label',{},'drive',{},'controller',{},'considerLosses',{},'pressure_rails',{}, ...
    'capArea',{},'rodArea',{},'ehaFixedDisplacement',{},'needsPressureOpt',{},'mAstar',{}, ...
    'seaStateIdx',{},'Hs',{},'Tp',{},'probability',{});
for fIdx = 1:size(familyDefs,1)
    [label, drive, controller, considerLosses, pressure_rails, needsPressureOpt, mAstarF] = familyDefs{fIdx,:};
    switch label
        case 'PassivePump', areaKey = 'PassivePump';
        case {'DHD2','DHD3','DHD4'}, areaKey = label;
        case {'EHA_mech','EHA_elec'}, areaKey = 'EHA';
    end
    for s = 1:nSeaStates
        t = struct('label',label,'drive',drive,'controller',controller,'considerLosses',considerLosses, ...
            'pressure_rails',pressure_rails,'capArea',sizedAreas.(areaKey).capArea, ...
            'rodArea',sizedAreas.(areaKey).rodArea,'ehaFixedDisplacement',strcmp(drive,'EHA'), ...
            'needsPressureOpt',needsPressureOpt,'mAstar',mAstarF,'seaStateIdx',s, ...
            'Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s),'probability',seaStates.probability(s));
        tasks(end+1) = t; %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d tasks (%d families x %d sea states).\n', nTasks, size(familyDefs,1), nSeaStates);

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
    outFile = fullfile(resultsDir, sprintf('task_%03d.mat', i));
    if isfile(outFile), continue; end

    runParams = struct('drive',task.drive,'controller',task.controller, ...
        'capArea',task.capArea,'rodArea',task.rodArea,'controlDT',CONTROL_DT);
    if ~isnan(task.considerLosses), runParams.considerLosses = task.considerLosses; end
    if ~isnan(task.pressure_rails), runParams.pressure_rails = task.pressure_rails; end
    if task.ehaFixedDisplacement
        runParams.ehaFixedDisplacement = true;
        runParams.shaftInertia = 249.9;
    end
    if strcmp(task.drive,'DHD')
        runParams.astarIterMax = ASTAR_ITER_MAX;
        runParams.mAstar = task.mAstar;
    end

    seaState = struct('Hs',task.Hs,'Tp',task.Tp);
    out = struct('label',task.label,'seaStateIdx',task.seaStateIdx,'Hs',task.Hs,'Tp',task.Tp, ...
        'probability',task.probability,'mAstar',task.mAstar);

    if task.needsPressureOpt
        runParams.highPressure = 35e6;
        res = optimizePressure(runParams, seaState, pressureGrid, struct(), denseMap);
        out.bestPressure = res.bestPressure; out.mechRGP = res.bestMechRGP; out.elecRGP = res.bestElecRGP;
        out.aveMechPow = res.bestAveMechPow; out.aveElecPow = res.bestAveElecPow;
        out.nAstarCapHits = res.totalAstarCapHits; out.fullGrid = res;
    else
        params = getParameters(runParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = task.Hs; params.simu.peakPeriod = task.Tp;
        wave = generateExcitingTorque(params);
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);
        ev = evaluate(params,dyn,ctrl);
        out.bestPressure = NaN; out.mechRGP = ev.mechRGP; out.elecRGP = ev.elecRGP;
        out.aveMechPow = ev.aveMechPow; out.aveElecPow = ev.aveElecPow; out.nAstarCapHits = ev.nAstarCapHits;
    end

    parsave(outFile, out);
    fprintf('[%d/%d] %s, sea state %d (Hs=%.2f,Tp=%.2f): mechRGP=%.3f elecRGP=%.3f\n', ...
        i, nTasks, task.label, task.seaStateIdx, task.Hs, task.Tp, out.mechRGP, out.elecRGP);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateFinalGridSearchResults();
plotFinalPressureSweeps();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
