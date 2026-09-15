% Run_PassivePump_SS1.m
% Single diagnostic run: PassivePump at Humboldt sea state 1, using its
% own optimal rail pressure for SS1 (15MPa, from
% results/mAstarConvAndPressureOpt/pressureOptSummary.csv, BestPressure_MPa
% column). Plots flap speed, flap torque, and mechanical power over time.
%
% Note: PassivePump's electrical loss (getValveLoss.m's PassivePump
% branch) is only tracked as a single whole-run scalar (eval.aveLoss),
% not a time series -- unlike EHA's copperCoeff*T_Act^2 loss, there's no
% instantaneous electrical-loss/power signal to plot here. "Power" below
% is mechanical power (eval.mechPower = -dyn.u.*dyn.thetaDot), available
% instantaneously; electrical output is a flat 85%-derated average only
% (ev.aveElecPow).
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

SEA_STATE_IDX = 1;
BEST_PRESSURE = 5e6; % Pa -- PassivePump's SS1-optimal rail pressure

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.PassivePump.capArea;
rodArea = S.sizedAreas.PassivePump.rodArea;

seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);

runParams = struct('drive','PassivePump','controller','CoulombDamping','pressure_rails',2, ...
    'capArea',capArea,'rodArea',rodArea,'controlDT',0.1,'highPressure',BEST_PRESSURE);

params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

fprintf('PassivePump, sea state %d (Hs=%.2f, Tp=%.2f), highPressure=%.0fMPa: mechRGP=%.4f elecRGP=%.4f\n', ...
    SEA_STATE_IDX, Hs, Tp, BEST_PRESSURE/1e6, ev.mechRGP, ev.elecRGP);
fprintf('aveMechPow=%.1f kW, aveElecPow=%.1f kW, aveLoss=%.1f kW\n', ...
    ev.aveMechPow/1e3, ev.aveElecPow/1e3, ev.aveLoss/1e3);

thetaDeg = dyn.theta * 180/pi;
fprintf('Flap position: max=%.2f deg, min=%.2f deg, peak |theta|=%.2f deg\n', ...
    max(thetaDeg), min(thetaDeg), max(abs(thetaDeg)));

%% Plot: flap speed, flap torque, mechanical power over time
fig = figure('Position',[100 100 1000 800],'Visible','off');

subplot(3,1,1)
plot(dyn.t, dyn.thetaDot, 'LineWidth', 0.75); grid on
ylabel('Flap speed (rad/s)'); title('Flap speed');

subplot(3,1,2)
plot(dyn.t, dyn.u/1e6, 'LineWidth', 0.75); grid on
ylabel('Flap torque (MN\cdotm)'); title('Flap torque');

subplot(3,1,3)
plot(dyn.t, ev.mechPower/1e3, 'LineWidth', 0.75); hold on
yline(ev.aveMechPow/1e3, '--', 'Color', [0.850 0.325 0.098], 'LineWidth', 1);
grid on
ylabel('Mechanical power (kW)'); xlabel('Time (s)');
title(sprintf('Mechanical power (dashed = ave %.1fkW)', ev.aveMechPow/1e3));
legend('Instantaneous','Average (post-ramp)','Location','best');

sgtitle(sprintf('PassivePump, sea state %d, %.0fMPa', SEA_STATE_IDX, BEST_PRESSURE/1e6));

outFile = fullfile(repoRoot,'results',sprintf('PassivePump_SS%d_%.0fMPa_speedTorquePower.png',SEA_STATE_IDX,BEST_PRESSURE/1e6));
print(fig, outFile, '-dpng', '-r150', '-painters');
fprintf('Wrote %s\n', outFile);
