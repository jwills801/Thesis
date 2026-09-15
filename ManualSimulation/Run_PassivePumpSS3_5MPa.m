% Run_PassivePumpSS3_5MPa.m
% Adds the one missing grid point results/mainResults/ needs for the
% "fall back to 5MPa when the symmetric +/-10deg position target is
% unreachable" rule: PassivePump SS3's fine grid ([14 19 21 23 25 27 29
% 31 35]MPa, centered on its 25MPa prior) never included 5MPa, unlike
% SS1's (which clips down to 5 since its 15MPa center is closer to the
% floor). Runs the identical task body Run_MainResults.m uses for a
% PassivePump pressureSweep task, and saves in the same
% task_<family>_ss<N>_p<MPa>.mat format so aggregateMainResults.m and
% plotMainResultsPressureSweeps.m pick it up like any other point.
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off, run once to backfill this single point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

resultsDir = fullfile(repoRoot,'results','mainResults');
SEA_STATE_IDX = 3;
PRESSURE = 5e6;

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);

runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
    'capArea',S.sizedAreas.PassivePump.capArea,'rodArea',S.sizedAreas.PassivePump.rodArea, ...
    'controlDT',0.1,'highPressure',PRESSURE);

params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

rampInd = round(params.simu.rampTime/params.simu.dt);
thetaPostRamp = dyn.theta(rampInd:end);
maxThetaDeg = max(thetaPostRamp)*180/pi;
minThetaDeg = min(thetaPostRamp)*180/pi;

out = struct('family','PassivePump','seaStateIdx',SEA_STATE_IDX,'Hs',Hs,'Tp',Tp, ...
    'pressure',PRESSURE,'mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
    'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'nAstarCapHits',ev.nAstarCapHits, ...
    'maxThetaDeg',maxThetaDeg,'minThetaDeg',minThetaDeg);

outFile = fullfile(resultsDir, sprintf('task_PassivePump_ss%d_p%05.1f.mat', SEA_STATE_IDX, PRESSURE/1e6));
save(outFile,'out');
fprintf('Wrote %s: mechRGP=%.4f elecRGP=%.4f max=%.2fdeg min=%.2fdeg\n', ...
    outFile, ev.mechRGP, ev.elecRGP, maxThetaDeg, minThetaDeg);
