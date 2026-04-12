function out = MPC_DP(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

fineTimeInd = length(uInd_history)+1;

% switch on timehorizon time intervals
if mod(fineTimeInd,ctrl.horizonInd) == 1
    % Find Which control index are in
    coarseTimeInd = floor((fineTimeInd-1)/ctrl.horizonInd);

% Get torque options
theta = states(2);
U = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
nU = length(U);
    % define nodes
    node = initlizeNodes(nU,states,ctrl.m_Astar);

    for uInd = 1:nU
        u = U(uInd);
        
    % Compute cost
        % Energy in this block
        w = ctrl.w;
        b_exc = ctrl.b_exc;
        Q_local = ctrl.Q_local;
        E_mech = (w'*x+b_exc)*u + u^2*Q_local;
        
        % Switching Loss        
        E_sw = 0;

        % Future Cost
        P = ctrl.P;
        rT = ctrl.rT(coarseTimeInd,:);
        b = ctrl.b(coarseTimeInd);
        E_term = x'*P*x + rT*x + b;

    % Choose the best option
    end

else
    uInd = uInd_history(end);
end
% Output 
out.controlValue = U(uInd);
out.controlIndex = uInd;
end

function node = initlizeNodes(nU,x0,m)
    node(nU) = struct();
    for uInd = 1:nU
        node(uInd).xf = x0;
        node(uInd).cost = 0;

        node(uInd).history = NaN(m,1);
        node(uInd).history(1) = uInd;
    end
end