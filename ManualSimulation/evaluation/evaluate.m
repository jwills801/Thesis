% evaluate.m
% Post-processes a completed run: drivetrain-specific losses (getEHALoss/
% getMotorLoss/getValveLoss), mechanical/electrical average power, the
% DHD/PassivePump flat mech-to-elec derating, and mechRGP/elecRGP
% (relative to ctrl.optTraj.avePow). getMotorLoss only applies to DHD
% with controller MPC_Astar_cont (its continuous electric trim) -- plain
% MPC_Astar has no such trim to charge a loss against.
% Calls: evaluation/getEHALoss.m, getMotorLoss.m (MPC_Astar_cont only), getValveLoss.m
% Called by: main_WEC_Simulation.m, optimization/optimizePressure.m,
%   diagnostics/checkEnergyBalance.m, validatePhase2Subset.m
function eval = evaluate(params,dyn,ctrl)
    
switch params.runParams.drive
    case 'DHD'
    % Calculate switching losses, plus (MPC_Astar_cont only) the series
    % hydraulic-to-electric conversion loss (main(1).tex "...Series
    % Hydraulic to Electric Converter") for its genuine continuous
    % electric trim on top of the discrete rail choice. Plain MPC_Astar
    % has no such trim -- u_elec=dyn.u-u_hyd is exactly 0 there (verified
    % directly) -- but getMotorLoss.m's EHA.LossFunc(Q,deltaP=0) is NOT
    % ~0 for nonzero flow Q, so calling it unconditionally was charging a
    % large (~tens of kW) fake loss for a conversion that isn't actually
    % happening. Confirmed to swing DHD2's elecRGP from ~0.48 to ~0.69 at
    % one sea state once removed -- this was a real, previously-uncaught
    % bug, not a rounding-level correction.
        eval = getValveLoss(params,dyn);
        if strcmp(params.runParams.controller,'MPC_Astar_cont')
            eval = getMotorLoss(eval,params,dyn);
        end
    case 'PassivePump'
    % Valve (check-valve) throttling loss only -- PassivePump has no
    % electric generator behind the valve. getMotorLoss.m's u_elec would
    % just be coulombDamping.m's tanh-smoothing artifact, not a real
    % actuator; running that through the EHA loss function was inflating
    % aveLoss by ~100x (1.4kW -> 142kW, measured). The main motor's flat
    % 90% conversion efficiency below still applies.
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
        % Main-motor efficiency 90% (was 85%) -- DHD/PassivePump's main
        % motor runs at steady, buffered conditions (unlike EHA's
        % continuously-varying operating point), so 85% was likely
        % conservative; see diagnostics/ReadMe.md.
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

    % Number of MPC_Astar control windows (0 for every other controller)
    % whose search hit its astarIterMax cap -- see control/MPC_Astar.m.
    % Surfaced here (rather than only via its console warning) so a whole
    % simulation's total is visible in the saved result, not just the log.
    eval.nAstarCapHits = dyn.nAstarCapHits;

end