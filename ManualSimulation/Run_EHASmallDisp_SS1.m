% Run_EHASmallDisp_SS1.m
% One-off experiment: what if the EHA's pump/motor unit displacement were
% 10x smaller (10.7 cc/rev instead of the real Pourmovahed et al. 1992b
% 107cc/rev unit both models/makeEHALossMap.m and
% makeEHALossMap_fixedDisp.m are built around)? Motivated by
% Run_EHAVarLossBreakdown_SS1.m showing EHA_var_elec's sea-state-1 loss is
% almost entirely a constant-speed "idle" friction loss
% (Scale*d*Cv*mu*w, independent of fracDisp) rather than copper loss.
%
% IMPORTANT caveat, confirmed algebraically before running this: Scale
% (how many 107cc/rev-equivalent units are ganged in parallel to hit the
% required max flow at 2000RPM) is defined as
% Scale = (maxFlow/w0)*2*pi*1e6/D*1.2, and d = D*1e-6/(2*pi) -- so
% Scale*d = 1.2*maxFlow/w0 EXACTLY, independent of D. Since the dominant
% idle-loss term (Scale*d*Cv*mu*w) depends on D only through that product,
% shrinking D while keeping the same total maxFlow (i.e. using
% proportionally more, smaller units -- the only physically consistent
% reading of "10x smaller displacement" here) leaves that dominant term
% UNCHANGED. Only the windage term (Scale*d^(5/3)*fracDisp*Ch*w^2/2,
% which scales as D^(2/3)) actually shrinks -- and it's already small
% here since fracDisp is tiny. Run anyway, for real numbers rather than
% just the algebra.
%
% Builds local D-parametrized copies of both EHA loss-map functions
% (D=10.7 instead of 107) and injects the result into params.hyd.EHA
% right after getParameters (same override pattern
% Run_MAstarConvAndPressureOpt.m uses for params.hyd.switchMap) -- does
% NOT modify models/makeEHALossMap*.m, so this is scoped to this script
% only.
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

SEA_STATE_IDX = 1;
D_SMALL = 10.7; % cc/rev, 10x smaller than the real 107cc/rev reference unit
copperCoeff = 5.0e-4;

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

seaStates = humboldtSeaStates();
Hs = seaStates.Hs(SEA_STATE_IDX); Tp = seaStates.Tp(SEA_STATE_IDX);
vMax = 1; % must match parameters/getHydraulic.m

cases = {
    'EHA_var_elec',   false
    'EHA_fixed_elec', true
    };

results = struct();
for c = 1:size(cases,1)
    label = cases{c,1}; fixedDisp = cases{c,2};

    runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',1, ...
        'capArea',capArea,'rodArea',rodArea,'controlDT',0.1,'ehaFixedDisplacement',fixedDisp);
    if fixedDisp
        runParams.shaftInertia = 249.9;
    end

    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

    % Override with the D=10.7cc/rev-based loss map (baseline uses
    % params.hyd.capArea, computed inside getParameters, as A)
    if fixedDisp
        params.hyd.EHA = makeEHALossMap_fixedDisp_D(params.hyd.capArea, vMax, copperCoeff, D_SMALL);
    else
        params.hyd.EHA = makeEHALossMap_D(params.hyd.capArea, vMax, copperCoeff, D_SMALL);
    end

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    fprintf('[D=%.1fcc/rev] %s, sea state %d: mechRGP=%.4f elecRGP=%.4f aveMechPow=%.1fkW aveElecPow=%.1fkW aveLoss=%.1fkW\n', ...
        D_SMALL, label, SEA_STATE_IDX, ev.mechRGP, ev.elecRGP, ev.aveMechPow/1e3, ev.aveElecPow/1e3, ev.aveLoss/1e3);

    results.(label) = struct('mechRGP',ev.mechRGP,'elecRGP',ev.elecRGP, ...
        'aveMechPow',ev.aveMechPow,'aveElecPow',ev.aveElecPow,'aveLoss',ev.aveLoss);

    if ~fixedDisp
        % Loss breakdown plot (hydraulic vs. electric/copper), same style
        % as Run_EHAVarLossBreakdown_SS1.m, but at D=10.7cc/rev
        [cap,~] = params.hyd.getVolandFlow(params,dyn.states);
        Q = cap.velA;
        theta = dyn.states(2,:)'; thetaDot = dyn.states(1,:)';
        F = dyn.u ./ params.hyd.Force2Torque(theta);
        deltaP = F/params.hyd.capArea;

        n = length(dyn.t);
        hydLoss = NaN(n,1); elecLoss = NaN(n,1); totLoss = NaN(n,1); fracDisp = NaN(n,1);
        for i = 1:n
            [totLoss(i), elecLoss(i), fracDisp(i)] = ehaVarLossSplit(Q(i),deltaP(i),vMax*capArea,copperCoeff,D_SMALL);
        end
        hydLoss = totLoss - elecLoss;
        mechPower = ev.mechPower;

        fig = figure('Position',[100 100 1000 700],'Visible','off');
        subplot(2,1,1)
        plot(dyn.t, hydLoss/1e3, 'LineWidth', 0.75); hold on
        plot(dyn.t, elecLoss/1e3, 'LineWidth', 0.75); grid on
        ylabel('Loss (kW)'); legend('Hydraulic (friction/windage/leakage)','Electric (copper)','Location','best');
        title(sprintf('EHA\\_var\\_elec loss breakdown, sea state %d, D=%.1fcc/rev (was 107)', SEA_STATE_IDX, D_SMALL));

        subplot(2,1,2)
        plot(dyn.t, mechPower/1e3, 'LineWidth', 0.75); hold on
        plot(dyn.t, mechPower/1e3 - totLoss/1e3, 'LineWidth', 0.75);
        yline(0,'k-');
        yline(ev.aveMechPow/1e3, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1);
        yline(ev.aveElecPow/1e3, '--', 'Color', [0.850 0.325 0.098], 'LineWidth', 1);
        grid on
        xlabel('Time (s)'); ylabel('Power (kW)');
        legend('Mechanical','Electrical (net)','','Ave mech','Ave elec','Location','best');
        title(sprintf('Power (ave elec = %.1fkW)', ev.aveElecPow/1e3));

        outFile = fullfile(repoRoot,'results',sprintf('EHA_var_elec_SS%d_lossBreakdown_D%.1fcc.png',SEA_STATE_IDX,D_SMALL));
        print(fig, outFile, '-dpng', '-r150', '-painters');
        fprintf('Wrote %s\n', outFile);
    end
end

fprintf('\n--- Before (D=107cc/rev, from results/EHA_fixedVsVariable/summary.csv, SS%d) vs after (D=%.1fcc/rev) ---\n', SEA_STATE_IDX, D_SMALL);
fprintf('EHA_var_elec:   before elecRGP=-0.6673 aveElecPow=-31.0kW  |  after elecRGP=%.4f aveElecPow=%.1fkW\n', ...
    results.EHA_var_elec.elecRGP, results.EHA_var_elec.aveElecPow/1e3);
fprintf('EHA_fixed_elec: before elecRGP=0.5147 aveElecPow=23.9kW  |  after elecRGP=%.4f aveElecPow=%.1fkW\n', ...
    results.EHA_fixed_elec.elecRGP, results.EHA_fixed_elec.aveElecPow/1e3);

%% Local D-parametrized copies of the loss-map builders
function EHA = makeEHALossMap_D(A,vMax,copperCoeff,D)
nw = 100; w_vals = linspace(-.5,.5,nw);
nT = 100; T_vals = linspace(-1,1,nT)*3e7;
[W,T] = ndgrid(w_vals,T_vals);
Loss = NaN(size(W));
for i = 1:length(W(:))
    V = W(i)*2.7574; Q = V*A;
    F = T(i)/2.7574; deltaP = F/A;
    Loss(i) = ehaLoss_D(Q,deltaP,vMax*A,copperCoeff,D);
end
w = reshape(W,[nw*nT,1]); t = reshape(T,[nw*nT,1])/1e6; l = reshape(Loss,[nw*nT,1])/1e5;
indep = [w.^2, w.*t, t.^2, ones(nT*nw,1)];
par = pinv(indep)*l;
coeffs = par*1e5;
EHA = struct();
EHA.LossCoeffs = [coeffs(1) coeffs(2)/1e6 coeffs(3)/1e6/1e6 coeffs(4)];
EHA.LossFunc = @(Q,deltaP) ehaLoss_D(Q,deltaP,vMax*A,copperCoeff,D);
end

function Loss = ehaLoss_D(Q,deltaP,maxFlow,copperCoeff,D)
Wrpm = 2000; w = Wrpm.*(2*pi/60);
Scale = maxFlow/w*2*pi*1e6/D *1.2;
d = (D*100^-3)/(2*pi);
Cf =  53.7e-3; Ch = 53.6; Cv = 23.5e3; Cs = 4.26e-9; Cst = 0*1e-5;
mu=(32e-6)*870; B = 1.7e9; rho = 870;

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
end

function [Loss, elecLoss, fracDisp] = ehaVarLossSplit(Q,deltaP,maxFlow,copperCoeff,D)
Wrpm = 2000; w = Wrpm.*(2*pi/60);
Scale = maxFlow/w*2*pi*1e6/D *1.2;
d = (D*100^-3)/(2*pi);
Cf =  53.7e-3; Ch = 53.6; Cv = 23.5e3; Cs = 4.26e-9; Cst = 0*1e-5;
mu=(32e-6)*870; B = 1.7e9; rho = 870;

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
elecLoss = P_L_elect;
end

function EHA = makeEHALossMap_fixedDisp_D(A,vMax,copperCoeff,D)
nw = 100; w_vals = linspace(-.5,.5,nw);
nT = 100; T_vals = linspace(-1,1,nT)*3e7;
[W,T] = ndgrid(w_vals,T_vals);
Loss = NaN(size(W));
for i = 1:length(W(:))
    V = W(i)*2.7574; Q = V*A;
    F = T(i)/2.7574; deltaP = F/A;
    Loss(i) = ehaLossFixed_D(Q,deltaP,vMax*A,copperCoeff,D);
end
w = reshape(W,[nw*nT,1]); t = reshape(T,[nw*nT,1])/1e6; l = reshape(Loss,[nw*nT,1])/1e5;
indep = [w.^2, w.*t, t.^2, ones(nT*nw,1)];
par = pinv(indep)*l;
coeffs = par*1e5;
EHA = struct();
EHA.LossCoeffs = [coeffs(1) coeffs(2)/1e6 coeffs(3)/1e6/1e6 coeffs(4)];
EHA.LossFunc = @(Q,deltaP) ehaLossFixed_D(Q,deltaP,vMax*A,copperCoeff,D);
end

function Loss = ehaLossFixed_D(Q,deltaP,maxFlow,copperCoeff,D)
Wrpm = 2000; w0 = Wrpm.*(2*pi/60);
Scale = maxFlow/w0*2*pi*1e6/D *1.2;
d = (D*100^-3)/(2*pi);
Cf =  53.7e-3; Ch = 53.6; Cv = 23.5e3; Cs = 4.26e-9; Cst = 0*1e-5;
mu=(32e-6)*870; B = 1.7e9; rho = 870;

fracDisp = 1;
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*w0*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
w = (Q + sign(deltaP)*QLoss) / (d*fracDisp*Scale);

TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );
T_Ideal = deltaP*d*fracDisp*Scale;
T_Act = T_Ideal + sign(w)*TLoss;

P_L_elect = copperCoeff*T_Act^2;
P_out = w*T_Act + P_L_elect;
P_in = Q*deltaP;
Loss = abs(P_in-P_out);
end
