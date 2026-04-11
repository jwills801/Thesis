function out = MPC_DP(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

timeInd = length(uInd_history)+1;

% Get position and velocity
thetaDot = states(1);
theta = states(2);

% Get torque options
ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);

% switch on timehorizon time intervals
if mod(timeInd,ctrl.horizonInd) == 1

    node(length(params.ptoOptions(:))) = struct();

    % Compute cost
        % Energy in this block
        E_mech = (w^T*x+b_exc)*u + u^2*Q;
        
        % Switching Loss        
        E_sw = 0;

        % Future Cost
        E_term = x^T*P*x + r^T*x + b;

    % Choose the best option

else
    uInd = uInd_history(end);
end
% Output 
out.controlValue = ptoTorqueOptions(uInd);
out.controlIndex = uInd;
end

