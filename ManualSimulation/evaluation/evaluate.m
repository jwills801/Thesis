function eval = evaluate(params,dyn)
    
switch params.runParams.drive
    case {'DHD', 'PassivePump'}
    % Calculate switching losses
        eval = getValveLoss(params,dyn);
        eval = getMotorLoss(eval,params,dyn);
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
        eval.aveElecPow = 0.85*eval.aveElecPow;
        eval.aveLoss = eval.aveLoss + 0.15*eval.aveElecPow;
    end

end