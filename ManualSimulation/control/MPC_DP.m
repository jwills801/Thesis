function out = MPC_DP(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control input indexes
if isempty(uInd_history), uIndPrev = 1; else uIndPrev = uInd_history(end); end

fineTimeInd = length(uInd_history)+1;

% Get torque options
theta = states(2);
nU = length(params.hyd.ptoForceOptions(:));
U = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);

% Only switch on timehorizon time intervals
if mod(fineTimeInd,ctrl.horizonInd) ~= 1
    uInd = uInd_history(end);
else
    % Find Which control index are in
    coarseTimeInd = floor((fineTimeInd-1)/ctrl.horizonInd);

    % define nodes
    node = initlizeNodes(nU,states);

    flag = 0; iter = 0; iterMax = 10;
    while flag == 0

    % sort nodes
    [~, idx] = sort([node.cost], 'descend');
    
    % reorder nodes
    node = node(idx);

    % Get info about this node
    k = length(node(1).history) + coarseTimeInd;
    uPrev = node(1).history(end);
    xf = node(1).xf;
    

    % look at all the possible options from this node
    costVec = NaN(nU,1); xfVec = NaN(4,nU);
    for uInd = 1:nU
        [costVec(uInd),xfVec(:,uInd)] = ComputeCost(ctrl,U(uInd),uPrev,xf,k);
    end
    % Pick best option
    [horizonCost,uInd] = max(costVec);

    % put info into this node
    node(1).cost = node(1).cost + horizonCost;
    node(1).history(end+1) = uInd;
    node(1).xf = xfVec(:,uInd);


        % Check if we reached the end of the horizon
            if k == ctrl.m_Astar
                flag = 1;
            end

    % Dont do too many iterations
    iter = iter +1;
    if iter>iterMax
        flag = -1;
    end
    end % while loop

end

% Output 
out.controlValue = U(uInd);
out.controlIndex = uInd;
end



function node = initlizeNodes(nU,x0)
    node(nU) = struct();
    for uInd = 1:nU
        node(uInd).xf = x0;
        node(uInd).cost = 0;
        node(uInd).history(1) = uInd;
    end
end


function [cost,xf] = ComputeCost(ctrl,u,uPrev,x,k)
    % Compute cost
        % Energy in this block
        w = ctrl.w;
        b_exc = ctrl.b_exc(k);
        Q_local = ctrl.Q_local;
        E_mech = (w'*x+b_exc)*u + u^2*Q_local;

        % Switching Loss
        E_sw = 0;

        % Future Cost
        P = ctrl.P;
        rT = ctrl.rT(k,:);
        b = ctrl.b(k);
        E_term = x'*P*x + rT*x + b;
        
        cost = E_mech + E_sw + E_term;

        % Advance state
        M_block = ctrl.M_block;
        T_block = ctrl.T_block(:,k);
        H_block = ctrl.H_block;
        xf = M_block*x+T_block+H_block*u;
end