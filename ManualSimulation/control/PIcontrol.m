function out = PIcontrol(params,ctrl,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

% Calculate the amount of time since the last switch
timeInd = length(uInd_history)+1;
lastSwitchInd = find(diff(uInd_history)~=0,1,'last');
timeSinceSwitch = params.simu.time(timeInd) - params.simu.time(lastSwitchInd);

% Get position and velocity
thetaDot = states(1);
theta = states(2);

% Get torque options
ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
if ctrl.limitChoices
    ptoTorqueOptions([5,6,8,9]) = NaN;
end

% if it hasnt been very long, use the previous control index
if timeSinceSwitch < .21
    uInd = uInd_history(end);
else
% if its been long enough, recalculate the control input
    H = freqresp(params.phys.sys , 2*pi/params.simu.peakPeriod);

    Kp = real(1/H(1)');
    Ki = .8*(-2*pi/5*imag(1/H(1)'));

    u_cont = -1*(Kp*thetaDot + Ki*theta);

    % Discretize  
    [~,uInd] = min(abs(u_cont-ptoTorqueOptions));
end

out.controlValue = ptoTorqueOptions(uInd);
out.controlIndex = uInd;
end