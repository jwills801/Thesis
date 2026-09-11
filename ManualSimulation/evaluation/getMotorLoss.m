% getMotorLoss.m
% DHD-only: models the "series hydraulic-to-electric converter" trim
% (the gap between the actual control torque and the nearest discrete
% rail torque, u_elec=dyn.u-u_hyd) as if it passed through an EHA-style
% generator, using params.hyd.EHA.LossFunc. Not called for PassivePump
% (that gap there is just coulombDamping.m's tanh-smoothing artifact, not
% a real actuator -- see evaluate.m's comment).
% Calls: none
% Called by: evaluation/evaluate.m
function eval = getMotorLoss(eval,params,dyn)

% Seperate hydraulic and electric torques
    % Hydraulic torques come from the valve operation
    % Electric torque comes from a series hydraulic to electric conversion
u_hyd = NaN(size(dyn.u));
for tInd = 1:length(dyn.t)
    % Hydraulic torque
    u_hyd(tInd) = params.hyd.Force2Torque(dyn.theta(tInd))*params.hyd.ptoForceOptions(dyn.uInd(tInd));
end
u_elec = dyn.u-u_hyd;

% Loss in the series hydraulic to electric conversion step
[~,rod] = params.hyd.getVolandFlow(params,dyn.states);
Q = rod.velA;
deltaP = u_elec ./ params.hyd.Force2Torque(dyn.states(2,:)')/ params.hyd.rodArea;

ElecPowloss = NaN(length(dyn.t),1);
for t_ind = 1:length(dyn.t)
    ElecPowloss(t_ind) = params.hyd.EHA.LossFunc(Q(t_ind),deltaP(t_ind));
end
ElecPowLossAfterRamp = ElecPowloss(dyn.t > params.simu.rampTime);

% Calculate Hydrualic and Electric powers
    % There should be no electric power if the h2e conversion is in series
powHyd = u_hyd.*dyn.thetaDot;
powElec = u_elec.*dyn.thetaDot;
powTot = dyn.u.*dyn.thetaDot;

% Cumulative energy 
E_hyd = cumtrapz(dyn.t,powHyd);
E_elec = cumtrapz(dyn.t,powElec);
E_tot = cumtrapz(dyn.t,powTot);

% Energy Loss from main motor to deal with ending acc vol
    % This doesnt use the EHA loss map because the size would be wrong
eff = .8;
Eloss_mm = (1-eff)*(-E_hyd(end)) * (E_hyd(end)<0);

% Total losses
eval.TotalLoss = eval.TotalLoss + trapz(dyn.t,ElecPowloss) + Eloss_mm;
eval.TotalLossAfterRamp = eval.TotalLossAfterRamp + sum(ElecPowLossAfterRamp)*params.simu.dt + Eloss_mm;
eval.aveLoss = eval.TotalLossAfterRamp / (params.simu.finalTime - params.simu.rampTime);

% Plots
if ~isfield(params.simu,'makePlots') || params.simu.makePlots
    figure, plot(dyn.t,E_hyd,dyn.t,E_elec,dyn.t,E_tot)
    legend('Hydraulic','Electric','Total'), ylabel('Energy [J]')
end
end
