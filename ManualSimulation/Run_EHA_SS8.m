% Run_EHA_SS8.m
% Single diagnostic run: EHA (fixed-displacement/variable-speed pump-motor,
% 20in bore from results/sizedAreas.mat) at a chosen Humboldt sea state and
% control mode (mechanical-only vs. electrical-loss-aware). Reports net
% electrical efficiency and plots (1) instantaneous efficiency over time,
% (2) a 4x1 subplot of flap torque, flap speed, hydraulic loss, and
% electric (copper) loss over time, and (3) mechanical vs. electrical
% power over time.
%
% Set SEA_STATE_IDX and CONSIDER_LOSSES below, then run. Output filenames
% are tagged with both (e.g. EHA_SS3_mech_power.png).
%
% Hydraulic/electric loss split: electricLoss(t) = copperCoeff*T_Act(t)^2
% (the actual I^2R copper heat, always >=0); hydraulicLoss(t) is the
% residual so the two sum exactly to the officially-tracked total loss
% (eval.loss, from models/makeEHALossMap_fixedDisp.m's ehaLoss). Both are
% recomputed here from the real (Q,deltaP,thetaDot) trajectory using the
% same physics/coefficients as that file.
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (top-level entry point)

SEA_STATE_IDX = 3;
CONSIDER_LOSSES = 0; % 1 = electrical-loss-aware (EHA_elec), 0 = mechanical-only (EHA_mech)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);
modeTag = 'elec'; if ~CONSIDER_LOSSES, modeTag = 'mech'; end
tag = sprintf('SS%d_%s', SEA_STATE_IDX, modeTag);

runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',CONSIDER_LOSSES, ...
    'capArea',capArea,'rodArea',rodArea,'controlDT',0.1, ...
    'ehaFixedDisplacement',true,'shaftInertia',249.9);

params = getParameters(runParams);
params.simu.makePlots = false;
params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);
dyn = timeLoop(params,wave,ctrl);
ev = evaluate(params,dyn,ctrl);

netEff = ev.aveElecPow/ev.aveMechPow;
fprintf('EHA_%s, sea state %d (Hs=%.2f, Tp=%.2f): mechRGP=%.4f elecRGP=%.4f\n', modeTag, SEA_STATE_IDX, Hs, Tp, ev.mechRGP, ev.elecRGP);
fprintf('aveMechPow=%.1f kW, aveElecPow=%.1f kW, aveLoss=%.1f kW, net efficiency=%.2f%%\n', ...
    ev.aveMechPow/1e3, ev.aveElecPow/1e3, ev.aveLoss/1e3, 100*netEff);

%% Recompute the hydraulic/electric loss split from the real trajectory
copperCoeff = 5.0e-4; % must match models/makeEHALossMap_fixedDisp.m's default
[cap,~] = params.hyd.getVolandFlow(params,dyn.states);
Q = cap.velA;
theta = dyn.states(2,:)';
thetaDot = dyn.states(1,:)';
F = dyn.u ./ params.hyd.Force2Torque(theta);
deltaP = F/params.hyd.capArea;

n = length(dyn.t);
hydLoss = NaN(n,1); elecLoss = NaN(n,1); totLoss = NaN(n,1);
for i = 1:n
    [totLoss(i), elecLoss(i)] = ehaLossSplit(Q(i),deltaP(i),1*capArea,copperCoeff);
end
hydLoss = totLoss - elecLoss;

mechPower = ev.mechPower;
instElecPower = mechPower - totLoss;
instEff = instElecPower ./ mechPower;

%% Plot 1: instantaneous efficiency over time
fig1 = figure('Position',[100 100 1000 500],'Visible','off');
plot(dyn.t, instEff, 'LineWidth', 0.75);
ylim([-2 1.2]);
grid on
xlabel('Time (s)'); ylabel('Instantaneous electrical efficiency');
title(sprintf('EHA %s instantaneous efficiency, sea state %d -- net efficiency = %.1f%%', modeTag, SEA_STATE_IDX, 100*netEff));

%% Plot 2: 4x1 subplot -- flap torque, flap speed, hydraulic loss, electric loss
fig2 = figure('Position',[100 100 1000 900],'Visible','off');
subplot(4,1,1)
plot(dyn.t, dyn.u/1e6, 'LineWidth', 0.75); grid on
ylabel('Flap torque (MN\cdotm)'); title('Flap torque');

subplot(4,1,2)
plot(dyn.t, thetaDot, 'LineWidth', 0.75); grid on
ylabel('Flap speed (rad/s)'); title('Flap speed');

subplot(4,1,3)
plot(dyn.t, hydLoss/1e3, 'LineWidth', 0.75); grid on
ylabel('Hydraulic loss (kW)'); title('Hydraulic (pump/motor friction+leakage) loss');

subplot(4,1,4)
plot(dyn.t, elecLoss/1e3, 'LineWidth', 0.75); grid on
ylabel('Electric loss (kW)'); title('Electric (copper) loss');
xlabel('Time (s)');

%% Plot 3: mechanical vs. electrical power over time
fig3 = figure('Position',[100 100 1000 500],'Visible','off');
plot(dyn.t, mechPower/1e3, 'LineWidth', 0.75); hold on
plot(dyn.t, instElecPower/1e3, 'LineWidth', 0.75);
yline(ev.aveMechPow/1e3, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1);
yline(ev.aveElecPow/1e3, '--', 'Color', [0.850 0.325 0.098], 'LineWidth', 1);
grid on
xlabel('Time (s)'); ylabel('Power (kW)');
legend('Mechanical power','Electrical power','Ave mechanical','Ave electrical','Location','best');
title(sprintf('EHA %s power over time, sea state %d (positive = generating)', modeTag, SEA_STATE_IDX));

print(fig1, fullfile(repoRoot,'results',sprintf('EHA_%s_instEfficiency.png',tag)), '-dpng', '-r150', '-painters');
print(fig2, fullfile(repoRoot,'results',sprintf('EHA_%s_torqueSpeedLoss.png',tag)), '-dpng', '-r150', '-painters');
print(fig3, fullfile(repoRoot,'results',sprintf('EHA_%s_power.png',tag)), '-dpng', '-r150', '-painters');

function [Loss, elecLoss] = ehaLossSplit(Q,deltaP,maxFlow,copperCoeff)
Wrpm = 2000;
w0 = Wrpm.*(2*pi/60);
Scale = maxFlow/w0*2*pi*1e6/107 *1.2;
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

fracDisp = 1;
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*w0*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
w = (Q + sign(deltaP)*QLoss) / (d*fracDisp*Scale);

TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );
T_Ideal = deltaP*d*fracDisp*Scale;
T_Act = T_Ideal + sign(w)*TLoss;

P_L_elect = copperCoeff*T_Act^2*sign(T_Act);
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;
Loss = abs(P_in-P_out);
elecLoss = copperCoeff*T_Act^2; % always-positive I^2R heat, for the plot split
end
