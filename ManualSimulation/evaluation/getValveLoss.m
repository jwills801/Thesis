function eval = getValveLoss(params,dyn)

%% Calculate volume and flow in each side
[cap,rod] = params.hyd.getVolandFlow(params,dyn.states);

switch params.runParams.drive
    case 'DHD' % Digital Hydraulic Drive
        % Find switching events
        % These denote the start of the switch
        % The start and end of the simulation are also counted as events
        eventInds = [1; find(diff(dyn.uInd)~=0); length(dyn.t)-1];
        eventTimes = dyn.t(eventInds);

        switchRate = length(eventInds)/params.simu.finalTime;

        switchMap = params.hyd.switchMap;
        mapDT = switchMap.finalTime;

        % Initilize the energy loss vector
        % This denotes the energy lost between event times
        % Thus there is one less loss entry than there are events
        loss = NaN(length(eventInds)-1,1);
        %
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
        % losses after ramp up
        lossAfterRamp = loss(eventTimes(1:end-1) > params.simu.rampTime);

        % output 
        eval.switchTimes = eventTimes;
        eval.aveSwitchRate = switchRate;
    case 'PassivePump'
        valveConstant = 1*params.hyd.capArea/sqrt(2e6);

        % Passive Pump only has steady loss
        cap.steadyLoss = openvalveLoss(params,cap,valveConstant,1,length(dyn.t)-1);
        rod.steadyLoss = openvalveLoss(params,rod,valveConstant,1,length(dyn.t)-1);

        loss = cap.steadyLoss + rod.steadyLoss;
        lossAfterRamp = loss;
end




% Output results
eval.TotalLoss = sum(loss);
eval.loss = loss;
eval.TotalLossAfterRamp = sum(lossAfterRamp);
eval.aveLoss = eval.TotalLossAfterRamp / (params.simu.finalTime - params.simu.rampTime);

% Optional plots
if 0
figure, plot(eventTimes(1:end-1),cap.switchLoss, ...
    eventTimes(1:end-1),rod.switchLoss,...
    eventTimes(1:end-1),cap.steadyLoss,...
    eventTimes(1:end-1),rod.steadyLoss)
legend('Cap Switch','Rod Switch','Cap Steady','Rod Steady')

figure, plot(eventTimes(1:end-1),cap.steadyLoss,...
    eventTimes(1:end-1),rod.steadyLoss)
legend('Cap Steady','Rod Steady')
end
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
