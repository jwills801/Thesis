function out = PIcontrol(params,ctrl,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

% Calculate the amount of time since the last switch
timeInd = length(uInd_history)+1;
lastSwitchInd = find(diff(uInd_history)~=0,1,'last');
timeSinceSwitch = params.simu.time(timeInd) - params.simu.time(lastSwitchInd);

% Get position and velocity
thetaDot = states(1);
theta = states(2);

% if its been long enough, recalculate the control input
H = freqresp(params.phys.sys , 2*pi/params.simu.peakPeriod);
Kp = real(1/H(1)');
Ki = .8*(-2*pi/5*imag(1/H(1)'));
u_cont = -1*(Kp*thetaDot + Ki*theta);

% if it hasnt been very long, and we are using DHD, then use the previous control index
switch params.runParams.drive
    case 'DHD'
        if timeSinceSwitch < .21
            uInd = uInd_history(end);
        else
            % Discretize
            ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
            [~,uInd] = min(abs(u_cont-ptoTorqueOptions));
        end
        out.controlValue = ptoTorqueOptions(uInd);
        out.controlIndex = uInd;
    case 'EHA'
        out.controlValue = u_cont;
        out.controlIndex = 1;
end


end