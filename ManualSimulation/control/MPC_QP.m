function out = MPC_QP(params,ctrl,wave,states,uInd_history)
% Get position and velocity
thetaDot = states(1);
theta = states(2);

% switch on timehorizon time intervals
timeInd = length(uInd_history)+1;
if mod(timeInd,ctrl.horizonInd) == 1
    % Horizon indices
    hInds = (timeInd:(timeInd+ctrl.MPC.m*ctrl.horizonInd-1));
    T = wave.torque.Texc(hInds);

    % X = M*states + H*T + H*L*u;
    % E = transpose(X)*C*u*.01;

    % Unpack matrices
    A = ctrl.MPC.A; Bx=ctrl.MPC.Bx; Bt = ctrl.MPC.Bt;
    x0 = states;
    B = x0'*Bx + T'*Bt;

    % Solve for u
    u = - (A+A') \ (B');

% If we are using EHA we use the continuous control, if DHD then we discretize
    switch params.runParams.drive
        case 'EHA'
            % for the eha we will use the control index to keep track of
            % our zero order hold
            out.controlValue = u(1);
            out.controlIndex = u(1);
        case 'DHD'
            % Get torque options
            ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
            [~,uInd] = min(abs(ptoTorqueOptions-u(1)));
            out.controlIndex = uInd;
            out.controlValue = ptoTorqueOptions(uInd);
    end
else % this is for if we are not recalculating the control, just using the previous value
        switch params.runParams.drive
        case 'EHA'
            % for the eha we will use the control index to keep track of
            % our zero order hold
            out.controlValue = uInd_history(end);
            out.controlIndex = uInd_history(end);
        case 'DHD'
            % Get torque options
            ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
            out.controlIndex = uInd_history(end);
            out.controlValue = ptoTorqueOptions(uInd_history(end));
        end
end
end

