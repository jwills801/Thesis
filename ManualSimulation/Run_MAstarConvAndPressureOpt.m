% Run_MAstarConvAndPressureOpt.m
% Combined supercomputer run, now that each family has its OWN locked
% bore diameter from the cylinder-area sweep (results/cylinderAreaSweep/
% summary.csv, sea state 8, 35MPa, mAstar=5): PassivePump=13in,
% DHD2=16in, DHD3=17in, DHD4=18in (see updateFamilyAreas.m, already run,
% results/sizedAreas.mat).
%
% Part A: A* search-horizon (m_Astar) convergence study at sea state 3,
%   for DHD2/DHD3/DHD4 (each at its own bore), fixed pressure -- but at
%   each family's OWN sea-state-3 near-optimal pressure (17/23/20 MPa for
%   DHD2/DHD3/DHD4 respectively, from Run_PressureOptSS3.m/
%   aggregatePressureOptSS3.m), not 35MPa (that's SS8's optimum -- a much
%   more energetic sea state -- not SS3's). Pressure is still fixed
%   across all mAstar values within a family, so this isolates the
%   horizon-length effect alone, same convention as
%   Run_AstarConvergenceStudy_100ms.m (which used sea state 7, the old
%   shared bore, and 35MPa -- appropriate there since SS7 sits much
%   closer to SS8 in energy than SS3 does).
% Part B: Pressure optimization (6-point grid, 5-35MPa) at a FIXED
%   m_Astar=5 for DHD, across all 4 families (PassivePump, DHD2, DHD3,
%   DHD4) and all 8 Humboldt sea states.
%
% ONE shared switch-loss map, sized to the LARGEST of the 3 DHD bores
% (18in, DHD4's own) -- valid for DHD2 (16in) and DHD3 (17in) too, since
% makeSwitchLossMap.m's actual loss physics depends only on real flow/
% volume at each grid point, never on capArea directly (capArea only
% sets how far the velA/vol axes extend) -- confirmed numerically in
% Run_CylinderAreaSweep.m's header/verifyMapIndependence.m check
% (0.41% agreement between two differently-sized maps at a shared query
% point). PassivePump needs no such map (its check-valve loss is sized
% directly from capArea/vMax inside getHydraulic.m).
%
% Task count: Part A = 3 families x 10 mAstar values = 30 single-sim
% tasks. Part B = 4 families x 8 sea states = 32 top-level tasks, each
% internally running a 6-point pressure grid via optimizePressure.m = 192
% underlying simulations. Total: 62 top-level (parfor) tasks, 222
% underlying simulation runs.
%
% Writes results/mAstarConvAndPressureOpt/task_astarConv_<family>_m<mAstar>.mat
% and task_pressureOpt_<family>_ss<seaStateIdx>.mat (one per task, saved
% immediately) and, via aggregateMAstarConvAndPressureOpt.m,
% results/mAstarConvAndPressureOpt/{astarConvSummary,pressureOptSummary}.csv.
%
% Calls: optimization/buildDenseSwitchMap.m, optimization/optimizePressure.m,
%   aggregateMAstarConvAndPressureOpt.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','mAstarConvAndPressureOpt');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;
MASTAR_FIXED = 5;
mAstarValues = [2 4 6 8 10 12 14 16 18 20];
pressureGrid = linspace(5e6, 35e6, 6);

dhdFamilies = {'DHD2',2; 'DHD3',3; 'DHD4',4};
allFamilies = {'PassivePump',NaN; 'DHD2',2; 'DHD3',3; 'DHD4',4};

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

%% Part A tasks: m_Astar convergence at sea state 3, each family at its
% OWN SS3-optimal pressure (from Run_PressureOptSS3.m/aggregatePressureOptSS3.m
% -- NOT 35MPa, which is SS8's optimum, not SS3's).
ss3Idx = 3;
Hs3 = seaStates.Hs(ss3Idx); Tp3 = seaStates.Tp(ss3Idx);
ss3Pressure = struct('DHD2',17e6, 'DHD3',23e6, 'DHD4',20e6);

tasks = struct('taskType',{},'family',{},'rails',{},'mAstar',{},'seaStateIdx',{},'pressure',{});
for f = 1:size(dhdFamilies,1)
    for m = mAstarValues
        tasks(end+1) = struct('taskType','astarConv','family',dhdFamilies{f,1}, ...
            'rails',dhdFamilies{f,2},'mAstar',m,'seaStateIdx',ss3Idx, ...
            'pressure',ss3Pressure.(dhdFamilies{f,1})); %#ok<AGROW>
    end
end
nAstarConvTasks = numel(tasks);

%% Part B tasks: pressure optimization at m_Astar=5, all families, all 8 sea states
% (pressure field unused here -- Part B finds its own best pressure per
% sea state via the internal grid search -- included only so this
% struct's fields match Part A's for concatenation)
for f = 1:size(allFamilies,1)
    for s = 1:nSeaStates
        tasks(end+1) = struct('taskType','pressureOpt','family',allFamilies{f,1}, ...
            'rails',allFamilies{f,2},'mAstar',MASTAR_FIXED,'seaStateIdx',s,'pressure',NaN); %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d top-level tasks (%d astarConv + %d pressureOpt, the latter each running a %d-point grid internally = %d underlying sims).\n', ...
    nTasks, nAstarConvTasks, nTasks-nAstarConvTasks, numel(pressureGrid), nAstarConvTasks + (nTasks-nAstarConvTasks)*numel(pressureGrid));

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

    if strcmp(task.taskType,'astarConv')
        outFile = fullfile(resultsDir, sprintf('task_astarConv_%s_m%02d.mat', task.family, task.mAstar));
        if isfile(outFile), continue; end

        runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
            'considerLosses',1,'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea, ...
            'highPressure',task.pressure,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',task.mAstar);
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

        out = struct('taskType','astarConv','family',task.family,'rails',task.rails,'mAstar',task.mAstar, ...
            'Hs',Hs3,'Tp',Tp3,'pressure',task.pressure,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
            'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits,'elapsedSec',elapsed);
        parsave(outFile, out);
        fprintf('[%d/%d] astarConv %s m=%d: mechRGP=%.3f elecRGP=%.3f capHits=%d (%.1f sec)\n', ...
            i, nTasks, task.family, task.mAstar, out.mechRGP, out.elecRGP, out.nAstarCapHits, elapsed);

    else % pressureOpt
        outFile = fullfile(resultsDir, sprintf('task_pressureOpt_%s_ss%d.mat', task.family, task.seaStateIdx));
        if isfile(outFile), continue; end

        Hs = seaStates.Hs(task.seaStateIdx); Tp = seaStates.Tp(task.seaStateIdx);
        probability = seaStates.probability(task.seaStateIdx);
        seaState = struct('Hs',Hs,'Tp',Tp);

        if strcmp(task.family,'PassivePump')
            runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
                'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea,'controlDT',CONTROL_DT);
            res = optimizePressure(runParams, seaState, pressureGrid, struct());
        else
            runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
                'considerLosses',1,'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea, ...
                'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',MASTAR_FIXED);
            res = optimizePressure(runParams, seaState, pressureGrid, struct(), sharedMap); %#ok<PFBNS>
        end

        out = struct('taskType','pressureOpt','family',task.family,'rails',task.rails,'mAstar',MASTAR_FIXED, ...
            'seaStateIdx',task.seaStateIdx,'Hs',Hs,'Tp',Tp,'probability',probability, ...
            'bestPressure',res.bestPressure,'mechRGP',res.bestMechRGP,'elecRGP',res.bestElecRGP, ...
            'aveMechPow',res.bestAveMechPow,'aveElecPow',res.bestAveElecPow,'nAstarCapHits',res.totalAstarCapHits, ...
            'fullGrid',res);
        parsave(outFile, out);
        fprintf('[%d/%d] pressureOpt %s ss%d: bestP=%.1fMPa mechRGP=%.3f elecRGP=%.3f capHits=%d\n', ...
            i, nTasks, task.family, task.seaStateIdx, out.bestPressure/1e6, out.mechRGP, out.elecRGP, out.nAstarCapHits);
    end
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateMAstarConvAndPressureOpt();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
