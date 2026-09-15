% Run_AstarConvergenceStudy_100ms.m
% A* search-horizon convergence study at the single most-annual-energy
% Humboldt sea state (sea state 7: Hs=3.66m, Tp=10.80s, 9.3% probability
% -- the highest probability-weighted theoretical power of the 8 Humboldt
% bins, beating both the highest-probability bin (3) and the largest-Hs
% bin (8); see the probability x avePow table computed alongside this
% script). Sweeps m_Astar (in 100ms coarse steps, so foresight time =
% 0.1*m_Astar) for DHD2, DHD3, and DHD4, all with the A* iteration cap
% raised to 100000 (per Run_DHD2_LargeCapStudy.m's fix -- confirmed to
% never actually fire at m_Astar=10 for DHD2). Pressure is fixed at
% 35 MPa for every task (not re-optimized per m_Astar) to isolate the
% horizon-length effect alone.
%
% mAstarValues goes up to 20 (2.0s foresight). DHD4's high-mAstar tasks in
% particular may be very slow or never converge in practice -- DHD4 has
% 16 force options per window (vs DHD2's 4), and the O(n^2) node-list
% resort in control/MPC_Astar.m scales badly with both branching factor
% and depth. That's fine given enough real cores (~30-32) to run every
% task in parallel: a straggler doesn't block the other 29, it just may
% not produce a file. Kill the session once progress has clearly stalled
% and run aggregateAstarConvergenceResults.m to recover everything that
% did finish.
%
% Requires results/denseSwitchMaps_100ms.mat (built by the DHD2 100ms
% probe in Run_DHD2_LargeCapStudy_100ms.m). ONE map covers all three
% families here -- DHD2/DHD3/DHD4 share identical capArea/rodArea
% (rail count is a control-layer choice, not a cylinder-sizing one; see
% results/sizedAreas.mat), and the map's PR axis is already a dense
% synthetic grid interpolated at runtime, not the actual few rails in
% use -- so the same map used for DHD2's 100ms study is exactly valid
% for DHD3/DHD4 too, no separate per-family build needed.
%
% Writes results/astarConvergence100ms/task_<family>_m<mAstar>.mat --
% ONE per task, saved immediately as it finishes (not only at the end of
% the parfor loop), so a task that completed before the session ends
% survives regardless. aggregateAstarConvergenceResults.m (separate,
% rerunnable) then builds summary.csv/.mat from whatever task files
% exist -- safe to call standalone at any time, including after an
% interrupted session, to recover partial results.
%
% Calls: aggregateAstarConvergenceResults.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','astarConvergence100ms');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

D = load(fullfile(repoRoot,'results','denseSwitchMaps_100ms.mat'));
denseMap = D.denseMapDHD2; % valid for all 3 families -- see header comment

% Fixed sea state: #7, the most annual-energy Humboldt bin (see header)
Hs = 3.65665236051502;
Tp = 10.8004291845494;
PRESSURE = 35e6;

CONTROL_DT = 0.1;
ASTAR_ITER_MAX = 100000;
mAstarValues = [2 4 6 8 10 12 14 16 18 20]; % 3 families x 10 values = 30 tasks -- with
    % ~30-32 real cores this runs fully in parallel (one worker per task), so a slow
    % DHD4-at-high-mAstar straggler doesn't delay the other 29; it just might not finish
    % (its per-task file simply won't exist -- see aggregateAstarConvergenceResults.m).
families = {'DHD2',2; 'DHD3',3; 'DHD4',4};

% Flatten into one task list (family x mAstar)
tasks = struct('family',{},'rails',{},'mAstar',{});
for f = 1:size(families,1)
    for m = mAstarValues
        tasks(end+1) = struct('family',families{f,1},'rails',families{f,2},'mAstar',m); %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('Built %d tasks (%d families x %d mAstar values) for sea state 7 (Hs=%.2f, Tp=%.2f).\n', ...
    nTasks, size(families,1), numel(mAstarValues), Hs, Tp);

nWorkers = min(nTasks, 4);
slurmCpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(slurmCpus)
    slurmCpus = str2double(getenv('SLURM_CPUS_ON_NODE'));
end
if ~isnan(slurmCpus) && slurmCpus > 0
    nWorkers = min(nTasks, slurmCpus);
end
try
    if isempty(gcp('nocreate'))
        parpool('local', nWorkers);
    end
catch ME
    fprintf('No parallel pool available (%s) -- running sequentially.\n', ME.message);
end

parfor i = 1:nTasks
    task = tasks(i); %#ok<PFBNS>
    outFile = fullfile(resultsDir, sprintf('task_%s_m%02d.mat', task.family, task.mAstar));
    if isfile(outFile)
        continue
    end
    runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',task.rails, ...
        'considerLosses',1,'capArea',sizedAreas.(task.family).capArea,'rodArea',sizedAreas.(task.family).rodArea, ...
        'highPressure',PRESSURE,'astarIterMax',ASTAR_ITER_MAX,'controlDT',CONTROL_DT,'mAstar',task.mAstar);
    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs;
    params.simu.peakPeriod = Tp;
    params.hyd.switchMap = denseMap;

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    tic;
    dyn = timeLoop(params,wave,ctrl);
    elapsed = toc;
    ev = evaluate(params,dyn,ctrl);

    out = struct('family',task.family,'rails',task.rails,'mAstar',task.mAstar,'Hs',Hs,'Tp',Tp, ...
        'pressure',PRESSURE,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP,'aveMechPow',ev.aveMechPow, ...
        'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits,'elapsedSec',elapsed);
    parsave(outFile, out);
    fprintf('[%d/%d] %s m=%d: mechRGP=%.3f elecRGP=%.3f capHits=%d (%.1f sec)\n', ...
        i, nTasks, task.family, task.mAstar, out.mechRGP, out.elecRGP, out.nAstarCapHits, elapsed);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateAstarConvergenceResults();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
