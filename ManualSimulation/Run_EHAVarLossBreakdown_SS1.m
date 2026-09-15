% Run_EHAVarLossBreakdown_SS1.m
% Loss breakdown for the worst variable-displacement EHA case from
% Run_EHA_FixedVsVariable.m: EHA_var_elec (loss-aware controller) at sea
% state 1, where results/EHA_fixedVsVariable/summary.csv shows
% elecRGP=-0.667 (aveMechPow=45.2kW but aveElecPow=-31.0kW -- net
% electrical CONSUMER, not generator).
%
% Same loss-split approach as Run_EHA_SS8.m (recompute hydraulic vs.
% electric/copper loss from the real (Q,deltaP,w) trajectory), but using
% the variable-displacement (constant 2000RPM shaft, chi solved per
% timestep) physics from models/makeEHALossMap.m's local ehaLoss,
% instead of the fixed-displacement one. Also plots fracDisp (chi) over
% time to show how small the commanded displacement fraction gets at
% this low-energy sea state, since models/makeEHALossMap.m's viscous
% loss term (Cv*mu*w) does NOT scale with fracDisp -- it's paid in full
% regardless of how little of the pump's displacement is actually being
% used, because the shaft never stops spinning at 2000RPM.
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

SEA_STATE_IDX = 1;

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);

runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',1, ...
    'capArea',capArea,'rodArea',rodArea,'controlDT',0.1,'ehaFixedDisplacement',false);

params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

fprintf('EHA_var_elec, sea state %d (Hs=%.2f, Tp=%.2f): mechRGP=%.4f elecRGP=%.4f\n', ...
    SEA_STATE_IDX, Hs, Tp, ev.mechRGP, ev.elecRGP);
fprintf('aveMechPow=%.1f kW, aveElecPow=%.1f kW, aveLoss=%.1f kW\n', ...
    ev.aveMechPow/1e3, ev.aveElecPow/1e3, ev.aveLoss/1e3);

%% Recompute the hydraulic/electric loss split (and fracDisp) from the real trajectory
copperCoeff = 5.0e-4; % must match models/makeEHALossMap.m's default
[cap,~] = params.hyd.getVolandFlow(params,dyn.states);
Q = cap.velA;
theta = dyn.states(2,:)';
thetaDot = dyn.states(1,:)';
F = dyn.u ./ params.hyd.Force2Torque(theta);
deltaP = F/params.hyd.capArea;

n = length(dyn.t);
hydLoss = NaN(n,1); elecLoss = NaN(n,1); totLoss = NaN(n,1); fracDisp = NaN(n,1);
for i = 1:n
    [totLoss(i), elecLoss(i), fracDisp(i)] = ehaLossSplit(Q(i),deltaP(i),1*capArea,copperCoeff);
end
hydLoss = totLoss - elecLoss;

mechPower = ev.mechPower;
instElecPower = mechPower - totLoss;

fprintf('Time-averaged (post-ramp) hyd loss=%.1fkW, elec(copper) loss=%.1fkW, |fracDisp| median=%.4f\n', ...
    mean(hydLoss(dyn.t>params.simu.rampTime))/1e3, mean(elecLoss(dyn.t>params.simu.rampTime))/1e3, ...
    median(abs(fracDisp(dyn.t>params.simu.rampTime))));

%% Plot: torque, speed, fracDisp, hyd vs elec loss, mech vs elec power
fig = figure('Position',[100 100 1100 1100],'Visible','off');

subplot(5,1,1)
plot(dyn.t, dyn.u/1e6, 'LineWidth', 0.75); grid on
ylabel('Flap torque (MN\cdotm)'); title('Flap torque');

subplot(5,1,2)
plot(dyn.t, thetaDot, 'LineWidth', 0.75); grid on
ylabel('Flap speed (rad/s)'); title('Flap speed');

subplot(5,1,3)
plot(dyn.t, fracDisp, 'LineWidth', 0.75); grid on
ylabel('fracDisp (\chi)'); title('Commanded displacement fraction (shaft always spins at 2000RPM regardless)');
ylim([-1.1 1.1])

subplot(5,1,4)
plot(dyn.t, hydLoss/1e3, 'LineWidth', 0.75); hold on
plot(dyn.t, elecLoss/1e3, 'LineWidth', 0.75); grid on
ylabel('Loss (kW)'); legend('Hydraulic (friction/windage/leakage)','Electric (copper)','Location','best');
title('Loss breakdown');

subplot(5,1,5)
plot(dyn.t, mechPower/1e3, 'LineWidth', 0.75); hold on
plot(dyn.t, instElecPower/1e3, 'LineWidth', 0.75);
yline(0,'k-');
yline(ev.aveMechPow/1e3, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1);
yline(ev.aveElecPow/1e3, '--', 'Color', [0.850 0.325 0.098], 'LineWidth', 1);
grid on
xlabel('Time (s)'); ylabel('Power (kW)');
legend('Mechanical','Electrical (net)','','Ave mech','Ave elec','Location','best');
title(sprintf('Power (dashed = averages; ave elec = %.1fkW, NEGATIVE = net consumer)', ev.aveElecPow/1e3));

sgtitle(sprintf('EHA\\_var\\_elec loss breakdown, sea state %d', SEA_STATE_IDX));

outFile = fullfile(repoRoot,'results',sprintf('EHA_var_elec_SS%d_lossBreakdown.png',SEA_STATE_IDX));
print(fig, outFile, '-dpng', '-r150', '-painters');
fprintf('Wrote %s\n', outFile);

function [Loss, elecLoss, fracDisp] = ehaLossSplit(Q,deltaP,maxFlow,copperCoeff)
% Variable-displacement (constant 2000RPM shaft) physics, copied from
% models/makeEHALossMap.m's local ehaLoss, split into hyd/elec (copper)
% components for plotting -- see that file for derivation/comments.
Wrpm = 2000;
w = Wrpm.*(2*pi/60);
Scale = maxFlow/w*2*pi*1e6/107 *1.2;
D = 107;
d = (D*100^-3)/(2*pi);
Cf =  53.7e-3;
Ch = 53.6;
Cv = 23.5e3;
Cs = 4.26e-9;
Cst = 0*1e-5;
mu=(32e-6)*870;
B = 1.7e9;
rho = 870;

fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale-sign(deltaP)*Scale*abs(d*w*deltaP/B));
if fracDisp <= 0
    fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale+sign(deltaP)*Scale*abs(d*w*deltaP/B));
end

T_Ideal = deltaP*d*fracDisp*Scale;
TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );
T_Act = T_Ideal + sign(w)*TLoss;

P_L_elect = copperCoeff*T_Act^2;
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;
Loss = abs(P_in-P_out);
elecLoss = P_L_elect; % always-positive I^2R heat, for the plot split
end
