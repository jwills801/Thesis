% evaluate.m
% Post-processes a completed run: drivetrain-specific losses (getEHALoss/
% getMotorLoss/getValveLoss), mechanical/electrical average power, the
% DHD/PassivePump flat mech-to-elec derating, and mechRGP/elecRGP
% (relative to ctrl.optTraj.avePow).
% Calls: evaluation/getEHALoss.m, getMotorLoss.m, getValveLoss.m
% Called by: main_WEC_Simulation.m, parameters/optimizePressure.m,
%   diagnostics/checkEnergyBalance.m, validatePhase2Subset.m
function eval = evaluate(params,dyn,ctrl)
    
switch params.runParams.drive
    case 'DHD'
    % Calculate switching losses, plus the series hydraulic-to-electric
    % conversion loss (main(1).tex "...Series Hydraulic to Electric
    % Converter") for controllers with a genuine continuous electric trim
    % (e.g. MPC_Astar_cont). For plain discrete controllers this trim is
    % ~0 so getMotorLoss contributes negligibly.
        eval = getValveLoss(params,dyn);
        eval = getMotorLoss(eval,params,dyn);
    case 'PassivePump'
    % Valve (check-valve) throttling loss only -- PassivePump has no
    % electric generator behind the valve. getMotorLoss.m's u_elec would
    % just be coulombDamping.m's tanh-smoothing artifact, not a real
    % actuator; running that through the EHA loss function was inflating
    % aveLoss by ~100x (1.4kW -> 142kW, measured). The main motor's flat
    % 85% conversion efficiency below still applies.
        eval = getValveLoss(params,dyn);
    case {'EHA'}
        eval = getEHALoss(params,dyn);
end

    % Calculate mechanical power
    eval.mechPower = -dyn.u.*dyn.thetaDot;

    % Calculate mechanical energy over time
    eval.mechEnergy = cumtrapz(dyn.t,eval.mechPower);

    % Average power after ramp period
    ind1 = params.simu.rampTime/params.simu.dt;
    ind2 = length(params.simu.time);
    inds = ind1:ind2;
    eval.aveMechPow = trapz(params.simu.time(inds),eval.mechPower(inds)) / (params.simu.finalTime - params.simu.rampTime);

    % Electrical output
    eval.aveElecPow = eval.aveMechPow - eval.aveLoss;

    switch params.runParams.drive
    case {'DHD', 'PassivePump'}
        % Compute the mech-to-elec conversion loss from the pre-derated
        % power, then derate: previously this used the already-derated
        % aveElecPow for the loss term, which silently broke the energy
        % balance (aveElecPow+aveLoss ~= aveMechPow) by ~2.25% of aveElecPow.
        mechToElecLoss = 0.15*eval.aveElecPow;
        eval.aveElecPow = eval.aveElecPow - mechToElecLoss;
        eval.aveLoss = eval.aveLoss + mechToElecLoss;
    end

    % Relative Generated Power, vs. the theoretical optimal average power
    % for this sea state (ctrl.optTraj.avePow). Moved here from plotAll.m
    % so RGP is available without having to call plotAll (which pops a
    % full figure set) -- needed for batch/sweep runs (Run_All_Cases.m,
    % pressure grid search).
    eval.mechRGP = eval.aveMechPow/ctrl.optTraj.avePow;
    eval.elecRGP = eval.aveElecPow/ctrl.optTraj.avePow;

end