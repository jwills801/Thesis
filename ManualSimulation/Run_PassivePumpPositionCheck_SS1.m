% Run_PassivePumpPositionCheck_SS1.m
% Independent re-verification of results/mainResults.csv's PassivePump
% SS1 row: reruns PassivePump at SS1's selected pressure (15MPa) and
% recomputes peak |theta| (post-ramp) using the EXACT same methodology as
% Run_MainResults.m (rampInd = round(rampTime/dt), theta taken from
% rampInd:end), to confirm the stored MaxThetaDeg/MinThetaDeg/
% MetPositionTarget=1 (peak|theta|=11.56deg, clears the 10deg target)
% weren't a fluke of that pipeline.
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off verification script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

SEA_STATE_IDX = 1;
SELECTED_PRESSURE = 15e6; % Pa -- results/mainResults.csv's PassivePump SS1 BestPressure_MPa

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.PassivePump.capArea;
rodArea = S.sizedAreas.PassivePump.rodArea;

seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);

runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
    'capArea',capArea,'rodArea',rodArea,'controlDT',0.1,'highPressure',SELECTED_PRESSURE);

params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

% Same methodology as Run_MainResults.m: post-ramp only
rampInd = round(params.simu.rampTime/params.simu.dt);
thetaPostRamp = dyn.theta(rampInd:end);
maxThetaDeg = max(thetaPostRamp)*180/pi;
minThetaDeg = min(thetaPostRamp)*180/pi;
peakThetaDeg = max(abs(maxThetaDeg), abs(minThetaDeg));

% Also report full-trajectory (including ramp-up) for comparison, since
% an earlier ad hoc check (Run_PassivePump_SS1.m) used the full
% trajectory rather than post-ramp -- want to see if that matters here.
thetaFullDeg = dyn.theta * 180/pi;
peakThetaFullDeg = max(abs(thetaFullDeg));

fprintf('PassivePump SS%d @ %.0fMPa: mechRGP=%.4f elecRGP=%.4f\n', SEA_STATE_IDX, SELECTED_PRESSURE/1e6, ev.mechRGP, ev.elecRGP);
fprintf('Post-ramp: max=%.2fdeg min=%.2fdeg peak|theta|=%.2fdeg\n', maxThetaDeg, minThetaDeg, peakThetaDeg);
fprintf('Full trajectory (incl. ramp-up): peak|theta|=%.2fdeg\n', peakThetaFullDeg);
fprintf('Clears 10deg target: %d   Clears 15deg target: %d\n', peakThetaDeg>10, peakThetaDeg>15);
fprintf('Cross-check against results/mainResults.csv: expected max=11.56 min=0.72 peak=11.56\n');
