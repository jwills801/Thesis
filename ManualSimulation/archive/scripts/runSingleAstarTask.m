function runSingleAstarTask(taskIdx)
% runSingleAstarTask.m
% Runs ONE row of results/astarHorizon/taskList.csv (built by
% prepareAstarTasks.m) and saves results/astarHorizon/task_NNN.mat.
% Meant to be invoked as its own OS process, e.g.
%   matlab -batch "runSingleAstarTask(7)"
% wrapped in `timeout` by run_overnight.sh, so a horizon length that's
% taking too long can simply be killed at the OS level -- this machine
% has no Parallel Computing Toolbox, so there's no way to cancel a
% running computation from inside a single MATLAB session. If a task is
% killed, it simply never writes its output file, which
% aggregateAstarHorizonResults.m treats the same as "not run yet".
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: run_overnight.sh (via matlab -batch, one process per task)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','astarHorizon');
outFile = fullfile(resultsDir, sprintf('task_%03d.mat', taskIdx));
if isfile(outFile)
    fprintf('Task %d already done, skipping.\n', taskIdx);
    return
end

taskList = readtable(fullfile(resultsDir,'taskList.csv'));
row = taskList(taskList.TaskIdx==taskIdx,:);
if height(row) ~= 1
    error('runSingleAstarTask:badTaskIdx','taskIdx %d not found (or duplicated) in taskList.csv',taskIdx);
end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.DHD2.capArea;
rodArea = S.sizedAreas.DHD2.rodArea;

D = load(fullfile(repoRoot,'results','denseSwitchMaps.mat'));
denseSwitchMap = D.denseSwitchMaps.DHD2;

runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',2, ...
    'considerLosses',1,'capArea',capArea,'rodArea',rodArea, ...
    'highPressure',row.HighPressure_Pa,'mAstar',row.MAstar);
params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = row.Hs;
params.simu.peakPeriod = row.Tp;
params.hyd.switchMap = denseSwitchMap;

fprintf('Task %d: sea state %d, m_Astar=%d, pressure=%.1fMPa -- starting\n', ...
    taskIdx, row.SeaStateIdx, row.MAstar, row.HighPressure_Pa/1e6);

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

out = struct('seaStateIdx',row.SeaStateIdx,'Hs',row.Hs,'Tp',row.Tp, ...
    'mAstar',row.MAstar,'highPressure',row.HighPressure_Pa, ...
    'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
    'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow);
save(outFile,'out');
fprintf('Task %d: DONE (mechRGP=%.3f elecRGP=%.3f)\n', taskIdx, ev.mechRGP, ev.elecRGP);
end
