function eval = getValveLoss(params,dyn,switchMap)

% Calculate volume and flow in each side
[cap,rod] = getVolandFlow(params,dyn);

% Find switching events
    % These denote the start of the switch
    % The start and end of the simulation are also counted as events
eventInds = [1, find(diff(dyn.uInd)~=0), length(dyn.t)-1];
eventTimes = dyn.t(eventInds);

switchRate = length(eventInds)/params.simu.finalTime;

mapDT = switchMap.finalTime;

% Initilize the energy loss vector
    % This denotes the energy lost between event times
    % Thus there is one less loss entry than there are events
loss = NaN(length(eventInds)-1,1);
%%
for k = 1:length(eventInds)-1
    % Which pressure rail did we switch from?
        % params.hyd.ptoTorqueOptions is a matrix
            % Each row is a different cap side option
            % Each col is a different rod side option
    [cap.switchFromInd,rod.switchFromInd] = ind2sub(size(params.hyd.ptoForceOptions),dyn.uInd(eventInds(k)));
    [cap.switchToInd,rod.switchToInd] = ind2sub(size(params.hyd.ptoForceOptions),dyn.uInd(eventInds(k+1)));

    % interpolate Losses
    cap.switchLoss(k) = interpolateLosses(params,switchMap,cap,eventInds(k));
    rod.switchLoss(k) = interpolateLosses(params,switchMap,rod,eventInds(k));

    % Penalize time between switches
    ind1 = eventInds(k) + round(mapDT/params.simu.dt)+1;
    ind2 = eventInds(k+1);

    % Find open valve loss from ind1 to ind2
    cap.steadyLoss(k) = openvalveLoss(params,cap,switchMap.valveConstant,ind1,ind2);
    rod.steadyLoss(k) = openvalveLoss(params,rod,switchMap.valveConstant,ind1,ind2);

end

% Sum up losses
    loss = cap.switchLoss + rod.switchLoss + cap.steadyLoss + rod.steadyLoss;

    figure, plot(eventTimes(1:end-1),cap.switchLoss, ...
        eventTimes(1:end-1),rod.switchLoss,...
        eventTimes(1:end-1),cap.steadyLoss,...
        eventTimes(1:end-1),rod.steadyLoss)
    legend('Cap Switch','Rod Switch','Cap Steady','Rod Steady')

        figure, plot(eventTimes(1:end-1),cap.steadyLoss,...
        eventTimes(1:end-1),rod.steadyLoss)
    legend('Cap Steady','Rod Steady')


% losses after ramp up
lossAfterRamp = loss(eventTimes(1:end-1) > params.simu.rampTime);


% Output results
eval.TotalValveLoss = sum(loss);
eval.lossAtSwitches = loss;
eval.switchTimes = eventTimes;
eval.aveSwitchRate = switchRate;
eval.TotalValveLossAfterRamp = sum(lossAfterRamp);
eval.aveValveLoss = eval.TotalValveLossAfterRamp / (params.simu.finalTime - params.simu.rampTime);
end

% Other functions called by getValveLoss
function [cap,rod] = getVolandFlow(params,dyn)
    theta = dyn.theta;
    thetaDot = dyn.thetaDot;
    
    % Calculate the length of the whole cylinder
    L = params.hyd.L(theta);
    
    % Calcuate the disctance from TDC
    xCap = L-params.hyd.L_retract;
        
    % Calulcate the distance from BDC
    xRod = params.hyd.stroke - xCap;

    % Calculate volumes in each side
    cap.vol = xCap*params.hyd.capArea;
    rod.vol = xRod*params.hyd.rodArea;

    % Calculate velocity of the cylinder
        % positive is in extension
    dLdt = params.hyd.dLdt(theta,thetaDot);

    % Calculate ideal flow into each cylinder
    cap.velA = dLdt*params.hyd.capArea;
    rod.velA = -dLdt*params.hyd.rodArea;
end

function switchLoss = interpolateLosses(params,switchMap,side,timeInd)
    % define variable to be interpolated on
    switchFrom = params.hyd.pressureRails(side.switchFromInd);
    switchTo = params.hyd.pressureRails(side.switchToInd);
    switchVelA = side.velA(timeInd);
    switchVol = side.vol(timeInd) + switchMap.hoseVolume;

    % interpolate
    switchLoss = interpn(switchMap.PR,switchMap.PR,switchMap.velA_vals,switchMap.vol_vals, switchMap.Eloss,...
        switchFrom,switchTo,switchVelA,switchVol);
end

function steadyLoss = openvalveLoss(params,side,k,timeInd1,timeInd2)
% Calculate the power loss in the valve while it is fully open
    powerLoss = (abs(side.velA(timeInd1:timeInd2))).^3/(k^2);

% Integrate to get the energy loss during this time
    if timeInd2>timeInd1
        steadyLoss = trapz(params.simu.time(timeInd1:timeInd2),powerLoss);
    else
        steadyLoss = 0;
    end
end
