function eval = evaluate(params,dyn)
    
    % load or calculate switching losses
    switchMap = makeSwitchLossMap(params);
        save("evaluation/SwitchMap.mat","switchMap")
    % load("SwitchMap.mat")

    % Calculate switching losses
    eval = getValveLoss(params,dyn,switchMap);

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
    eval.aveElecPow = eval.aveMechPow - eval.aveValveLoss;

end