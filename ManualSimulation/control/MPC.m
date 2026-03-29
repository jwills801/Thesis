function out = MPC(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

% Calculate the amount of time since the last switch
timeInd = length(uInd_history)+1;
lastSwitchInd = find(diff(uInd_history)~=0,1,'last');
    if isempty(lastSwitchInd)
        timeSinceSwitch = 10;
        uInd_history = 1;
    else
        timeSinceSwitch = params.simu.time(timeInd) - params.simu.time(lastSwitchInd);
    end

% Get position and velocity
thetaDot = states(1);
theta = states(2);

% Get torque options
ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);

% switch on timehorizon time intervals
if mod(timeInd,ctrl.horizonInd) == 1

% pad vectors so we can simulate past the official time vector
Texc = [wave.torque.Texc;zeros(ctrl.MPC.m*ctrl.horizonInd,1)];

% place holder u values
    % u = ptoTorqueOptions(ones(m,1));

    % Horizon indices
    hInds = (timeInd:(timeInd+ctrl.MPC.m*ctrl.horizonInd-1));

    % X = M*states + H*Texc(hInds) + H*L*u;

    % E = transpose(X)*C*u*.01;

    f_sw = zeros(4*ctrl.MPC.m*ctrl.horizonInd,1); f_sw(1) = -ctrl.MPC.gamma*2*ptoTorqueOptions(uInd_history(end));
    
    u = inv(-ctrl.MPC.Q+transpose(-ctrl.MPC.Q)) * transpose(ctrl.MPC.C)*(ctrl.MPC.M*states+ctrl.MPC.H*Texc(hInds)+f_sw);

    [~,uInd] = min(abs(ptoTorqueOptions-u(1)));
else
    uInd = uInd_history(end);
end
% Output 
out.controlValue = ptoTorqueOptions(uInd);
out.controlIndex = uInd;
end

