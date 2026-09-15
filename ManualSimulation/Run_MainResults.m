% Run_MainResults.m
% Combined supercomputer run intended to produce the thesis's main
% results, across all 8 Humboldt sea states:
%
%  - PassivePump: pressure-optimized per sea state, preferring a pressure
%    that lets flap position exceed 10deg (peak |theta|, post-ramp) when
%    the grid offers one, otherwise falling back to whichever pressure
%    maximizes elecRGP.
%  - DHD2/DHD3/DHD4: pressure-optimized per sea state (plain
%    argmax(elecRGP), no position constraint), mAstar=5 (matches the
%    prior sweep).
%  - EHA (fixed-displacement, electric-loss-aware only -- no mechanical-
%    only variant, per standing preference): single run per sea state,
%    no pressure grid (EHA has no pressure rails).
%
% Pressure grid, per (family, sea state): where results/mAstarConvAndPressureOpt/
% pressureOptSummary.csv already has an approximate optimum for that exact
% combo (plus DHD4 SS4=20MPa, confirmed separately -- not in that CSV),
% use a FINE grid centered on it: center + [-11,-6,-4,-2,0,2,4,6,11] MPa,
% clipped to [5,35]MPa and deduplicated. Where there's no prior at all --
% DHD2 SS4, DHD3 SS5-8, DHD4 SS2/3/5/6/7/8 -- those combos are missing
% because something about them makes the search take a very long time
% (not just unfinished bookkeeping), so they intentionally keep the OLD
% COARSE grid ([5 11 17 23 29 35]MPa, 6pts) rather than the finer one, to
% avoid multiplying the number of potentially-slow points.
%
% Hang robustness (isolation-only, no cancellation): every
% (family, sea state, pressure) triple -- not just (family, sea state) --
% is its own top-level parfor task, writing its own task_*.mat the moment
% it finishes (same if-isfile-continue resumability as every other run
% script here). If one point hangs, it only ties up the one worker
% running it for the rest of this job's walltime -- every other point,
% for every other combo, still finishes and lands on disk independently.
% [RUNNING]/[DONE] lines (with elapsed seconds) are printed around each
% task so a [RUNNING] line with no matching [DONE] in the log is exactly
% the hung-task list. Aggregation is intentionally a SEPARATE script
% (aggregateMainResults.m), not auto-chained here, so it can be run at
% any time against whatever task files exist -- including while this
% sweep is still running or stuck on a few hung points.
%
% Writes results/mainResults/task_<family>_ss<N>_p<MPa>.mat (PassivePump/
% DHD) and task_EHA_fixed_elec_ss<N>.mat (EHA). Run aggregateMainResults.m
% separately afterward (or anytime) to collect results/mainResults.csv.
%
% Calls: optimization/buildDenseSwitchMap.m
% Called by: none (top-level entry point, meant for cluster submission)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','mainResults');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;
MASTAR_FIXED = 5; % kept at 5, matching the prior sweep, per explicit confirmation

COARSE_GRID = [5 11 17 23 29 35] * 1e6; % Pa -- unchanged fallback for no-prior combos
FINE_OFFSETS = [-11 -6 -4 -2 0 2 4 6 11] * 1e6; % Pa, added to a known-prior center

POSITION_TARGET_DEG = 10; % PassivePump: prefer pressures where peak|theta| exceeds this

dhdFamilies = {'DHD2',2; 'DHD3',3; 'DHD4',4};

%% Load approximate-optimum priors from the earlier coarse sweep
priorFile = fullfile(repoRoot,'results','mAstarConvAndPressureOpt','pressureOptSummary.csv');
priorTable = readtable(priorFile);

priors = containers.Map('KeyType','char','ValueType','double');
for i = 1:height(priorTable)
    key = sprintf('%s_%d', priorTable.Family{i}, priorTable.SeaStateIdx(i));
    priors(key) = priorTable.BestPressure_MPa(i) * 1e6; % Pa
end
priors('DHD4_4') = 20e6; % confirmed separately, not in the CSV

%% Build (or load) the shared 100ms switch-loss map (DHD only), sized to DHD4's bore
mapFile = fullfile(repoRoot,'results','denseSwitchMap_shared_100ms.mat');
if isfile(mapFile)
    fprintf('Loading cached shared switch map from %s\n', mapFile);
    D = load(mapFile);
    sharedMap = D.denseMap;
else
    fprintf('Building 100ms dense switch map sized to DHD4''s bore (capArea=%.6f)...\n', sizedAreas.DHD4.capArea);
    tic;
    hydTmp = struct('capArea',sizedAreas.DHD4.capArea,'rodArea',sizedAreas.DHD4.rodArea,'stroke',5);
    denseMap = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.1); %#ok<NASGU>
    fprintf('Build took %.1f minutes\n', toc/60);
    save(mapFile,'denseMap');
    sharedMap = denseMap;
end

%% Build the flat task list (Hs/Tp resolved here, not inside parfor)
tasks = struct('type',{},'family',{},'rails',{},'seaStateIdx',{},'Hs',{},'Tp',{},'pressure',{});

% PassivePump + DHD2/DHD3/DHD4: one task per (family, sea state, pressure)
allRailFamilies = [{'PassivePump',NaN}; dhdFamilies];
for f = 1:size(allRailFamilies,1)
    family = allRailFamilies{f,1}; rails = allRailFamilies{f,2};
    for ss = 1:nSeaStates
        key = sprintf('%s_%d', family, ss);
        if isKey(priors,key)
            grid = unique(min(max(priors(key) + FINE_OFFSETS, 5e6), 35e6));
        else
            grid = COARSE_GRID;
        end
        for p = grid
            tasks(end+1) = struct('type','pressureSweep','family',family,'rails',rails, ...
                'seaStateIdx',ss,'Hs',seaStates.Hs(ss),'Tp',seaStates.Tp(ss),'pressure',p); %#ok<AGROW>
        end
    end
end

% EHA fixed-displacement, electric-loss-aware: one task per sea state, no grid
for ss = 1:nSeaStates
    tasks(end+1) = struct('type','ehaFixedElec','family','EHA_fixed_elec','rails',NaN, ...
        'seaStateIdx',ss,'Hs',seaStates.Hs(ss),'Tp',seaStates.Tp(ss),'pressure',NaN); %#ok<AGROW>
end

nTasks = numel(tasks);
fprintf('\nBuilt %d top-level tasks.\n', nTasks);

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
    Hs = task.Hs; Tp = task.Tp;

    if strcmp(task.type,'ehaFixedElec')
        label = sprintf('%s SS%d', task.family, task.seaStateIdx);
        outFile = fullfile(resultsDir, sprintf('task_EHA_fixed_elec_ss%d.mat', task.seaStateIdx));
        if isfile(outFile), continue; end
        fprintf('[RUNNING] %s\n', label);
        taskTimer = tic;

        runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',1, ...
            'capArea',sizedAreas.EHA.capArea,'rodArea',sizedAreas.EHA.rodArea, ...
            'controlDT',CONTROL_DT,'ehaFixedDisplacement',true,'shaftInertia',249.9);
        params = getParameters(runParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

        wave = generateExcitingTorque(params);
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);
        ev = evaluate(params,dyn,ctrl);

        out = struct('family',task.family,'seaStateIdx',task.seaStateIdx,'Hs',Hs,'Tp',Tp, ...
            'pressure',NaN,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
            'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'nAstarCapHits',0, ...
            'maxThetaDeg',NaN,'minThetaDeg',NaN);
        parsave(outFile, out);

        elapsed = toc(taskTimer);
        fprintf('[DONE]    %s (%.1fs) mechRGP=%.3f elecRGP=%.3f\n', label, elapsed, out.mechRGP, out.elecRGP);
        continue
    end

    % pressureSweep task: PassivePump or DHD2/3/4
    pMPa = task.pressure/1e6;
    label = sprintf('%s SS%d p=%.1fMPa', task.family, task.seaStateIdx, pMPa);
    outFile = fullfile(resultsDir, sprintf('task_%s_ss%d_p%05.1f.mat', task.family, task.seaStateIdx, pMPa));
    if isfile(outFile), continue; end
    fprintf('[RUNNING] %s\n', label);
    taskTimer = tic;

    if strcmp(task.family,'PassivePump')
        runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
            'capArea',sizedAreas.PassivePump.capArea,'rodArea',sizedAreas.PassivePump.rodArea, ...
            'controlDT',CONTROL_DT,'highPressure',task.pressure);
    else
        runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
            'considerLosses',1,'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea, ...
            'highPressure',task.pressure,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',MASTAR_FIXED);
    end

    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;
    if strcmp(task.family,'DHD2') || strcmp(task.family,'DHD3') || strcmp(task.family,'DHD4')
        params.hyd.switchMap = sharedMap; %#ok<PFBNS>
    end

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    rampInd = round(params.simu.rampTime/params.simu.dt);
    thetaPostRamp = dyn.theta(rampInd:end);
    maxThetaDeg = max(thetaPostRamp)*180/pi;
    minThetaDeg = min(thetaPostRamp)*180/pi;

    out = struct('family',task.family,'seaStateIdx',task.seaStateIdx,'Hs',Hs,'Tp',Tp, ...
        'pressure',task.pressure,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
        'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits, ...
        'maxThetaDeg',maxThetaDeg,'minThetaDeg',minThetaDeg);
    parsave(outFile, out);

    elapsed = toc(taskTimer);
    fprintf('[DONE]    %s (%.1fs) mechRGP=%.3f elecRGP=%.3f peak|theta|=%.1fdeg\n', ...
        label, elapsed, out.mechRGP, out.elecRGP, max(abs(maxThetaDeg),abs(minThetaDeg)));
end

fprintf('\nAll tasks submitted/completed. Run aggregateMainResults.m separately to collect results.\n');

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
