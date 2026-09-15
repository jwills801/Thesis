% Run_EHA_AllSeaStates.m
% EHA across all 8 Humboldt sea states, both control variants:
%   EHA_mech: MPC_QP optimized for mechanical power only (considerLosses=0)
%   EHA_elec: MPC_QP optimized for electrical power, loss-aware (considerLosses=1)
% Fixed-displacement/variable-speed pump-motor (chi=1, shaft speed solved),
% 20in bore from results/sizedAreas.mat, controlDT=0.1 (100ms). No pressure
% grid search needed -- EHA's control is continuous, not discrete-rail.
%
% Requires, already current as of this script's first run:
%   - results/sizedAreas.mat (EHA capArea/rodArea)
%
% Writes results/EHA_allSeaStates/task_%03d.mat (one per task, saved the
% moment it finishes) and, via aggregateEHAAllSeaStates.m,
% results/EHA_allSeaStates/summary.csv.
%
% Calls: wave/humboldtSeaStates.m, aggregateEHAAllSeaStates.m
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','EHA_allSeaStates');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

CONTROL_DT = 0.1;

% Sim length -- full length now (getSimulation.m's defaults). The short
% 100s/20s-ramp sanity pass earlier showed SS7/SS8 (longest periods)
% aren't converged at that length -- this run uses full length to get a
% trustworthy, current-parameter (shaftInertia=249.9, copperCoeff=5.0e-4)
% picture across all 8 sea states.
FINAL_TIME = 500; % s
RAMP_TIME = 50; % s

%% Build the flat task list: {label, considerLosses} x sea states
familyDefs = {
    'EHA_mech', 0
    'EHA_elec', 1
    };

tasks = struct('label',{},'considerLosses',{},'seaStateIdx',{},'Hs',{},'Tp',{},'probability',{});
for fIdx = 1:size(familyDefs,1)
    [label, considerLosses] = familyDefs{fIdx,:};
    for s = 1:nSeaStates
        t = struct('label',label,'considerLosses',considerLosses,'seaStateIdx',s, ...
            'Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s),'probability',seaStates.probability(s));
        tasks(end+1) = t; %#ok<AGROW>
    end
end
nTasks = numel(tasks);
fprintf('\nBuilt %d tasks (%d EHA variants x %d sea states).\n', nTasks, size(familyDefs,1), nSeaStates);

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

    runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',task.considerLosses, ...
        'capArea',capArea,'rodArea',rodArea,'controlDT',CONTROL_DT, ...
        'ehaFixedDisplacement',true,'shaftInertia',249.9);

    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = task.Hs; params.simu.peakPeriod = task.Tp;
    params.simu.finalTime = FINAL_TIME; params.simu.rampTime = RAMP_TIME;
    params.simu.time = (0:params.simu.dt:params.simu.finalTime)';
    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    out = struct('label',task.label,'seaStateIdx',task.seaStateIdx,'Hs',task.Hs,'Tp',task.Tp, ...
        'probability',task.probability,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
        'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'aveLoss',ev.aveLoss);

    parsave(outFile, out);
    fprintf('[%d/%d] %s, sea state %d (Hs=%.2f,Tp=%.2f): mechRGP=%.3f elecRGP=%.3f\n', ...
        i, nTasks, task.label, task.seaStateIdx, task.Hs, task.Tp, out.mechRGP, out.elecRGP);
end

fprintf('\nAll tasks submitted/completed. Aggregating...\n');
aggregateEHAAllSeaStates();

function parsave(outFile, out) %#ok<INUSD>
save(outFile,'out');
end
