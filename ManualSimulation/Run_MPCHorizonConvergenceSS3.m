% Run_MPCHorizonConvergenceSS3.m
% QP look-ahead-horizon convergence study for the EHA's electrical-loss-
% aware controller (MPC_QP, considerLosses=1) at sea state 3. Sweeps
% runParams.mpcHorizonPeriods (new override added to getControl.m's
% MPC_QP case -- default 2.5, the number of wave periods the receding-
% horizon QP looks ahead each time it re-solves) to check whether the
% current default is already converged, or could be shortened (cheaper)
% or needs lengthening.
%
% Uses EHA's own locked 20in bore, current shaftInertia=249.9/
% copperCoeff=5.0e-4. Full length (500s/50s ramp).
%
% Writes results/mpcHorizonConvSS3/task_h<horizonPeriods>.mat (one per
% task, saved immediately) and results/mpcHorizonConvSS3/summary.csv.
%
% Calls: none
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','mpcHorizonConvSS3');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

seaStates = humboldtSeaStates();
ss3Idx = 3;
Hs3 = seaStates.Hs(ss3Idx); Tp3 = seaStates.Tp(ss3Idx);

horizonPeriodsValues = [0.25 0.5 0.75 1 1.5 2 2.5 3 4 5 6 8];

fprintf('MPC_QP horizon convergence, EHA_elec, sea state 3 (Hs=%.2f, Tp=%.2f)\n', Hs3, Tp3);
fprintf('%-16s %-10s %-10s %-14s %-14s %-10s\n', 'horizonPeriods', 'mechRGP', 'elecRGP', 'aveMechPow_kW', 'aveElecPow_kW', 'elapsed_s');

rows = cell(numel(horizonPeriodsValues),1);
for i = 1:numel(horizonPeriodsValues)
    hp = horizonPeriodsValues(i);
    outFile = fullfile(resultsDir, sprintf('task_h%g.mat', hp));
    if isfile(outFile)
        d = load(outFile); out = d.out;
    else
        runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',1, ...
            'capArea',capArea,'rodArea',rodArea,'controlDT',0.1, ...
            'ehaFixedDisplacement',true,'shaftInertia',249.9,'mpcHorizonPeriods',hp);
        params = getParameters(runParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = Hs3;
        params.simu.peakPeriod = Tp3;

        wave = generateExcitingTorque(params);
        tic;
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);
        elapsed = toc;
        ev = evaluate(params,dyn,ctrl);

        out = struct('horizonPeriods',hp,'numHorizons',ctrl.numHorizons,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
            'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'elapsedSec',elapsed);
        save(outFile,'out');
    end
    fprintf('%-16g %-10.4f %-10.4f %-14.2f %-14.2f %-10.1f\n', ...
        out.horizonPeriods, out.mechRGP, out.elecRGP, out.aveMechPow/1e3, out.aveElecPow/1e3, out.elapsedSec);
    rows{i} = table(out.horizonPeriods, out.numHorizons, out.mechRGP, out.elecRGP, ...
        out.aveMechPow/1e3, out.aveElecPow/1e3, out.elapsedSec, ...
        'VariableNames', {'HorizonPeriods','NumHorizons','MechRGP','ElecRGP','AveMechPow_kW','AveElecPow_kW','ElapsedSec'});
end

summary = vertcat(rows{:});
writetable(summary, fullfile(resultsDir,'summary.csv'));
save(fullfile(resultsDir,'summary.mat'),'summary');
fprintf('\nWrote %s\n', fullfile(resultsDir,'summary.csv'));
