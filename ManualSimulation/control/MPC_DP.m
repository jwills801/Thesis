function out = MPC_DP(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control input indexes
if isempty(uInd_history), uIndPrev = 1; else uIndPrev = uInd_history(end); end

fineTimeInd = length(uInd_history)+1;

% Get torque options
theta = states(2);
nU = length(params.hyd.ptoForceOptions(:));

% Only switch on timehorizon time intervals
if mod(fineTimeInd,ctrl.horizonInd) ~= 1
    uInd = uInd_history(end);
else
    % Find Which control index are in
    k = floor((fineTimeInd-1)/ctrl.horizonInd)+1;

    % define nodes
    nodes = initlizeNodes(nU,states,params,ctrl,k,uIndPrev);

    % Begin Astar algorith 
    flag = 0; iter = 0; iterMax = 10;
    while flag == 0

        % sort nodes
        [~, idx] = sort([nodes.cost], 'ascend');
        
        % reorder nodes
        nodes = nodes(idx);

        % Get info about this node
        d = length(nodes(1).history);
        uPrev = nodes(1).history(end);
        xf = nodes(1).xf;

        % Check if we reached the end of the horizon
        if d == ctrl.m_Astar
            flag = 1;
            break
        end


        % look at all the possible options from this node
        for uInd = 1:nU
            [E_step, E_term, x_next] = ComputeCost(params, ctrl, uInd, uPrev, xf, k, d+1);
        
            new_node.history = [nodes(1).history, uInd];
            new_node.xf = x_next;
            new_node.E_absorbed = nodes(1).E_absorbed + E_step;
            new_node.cost = nodes(1).E_absorbed + E_term;

            % Add child to the pool
            nodes = [nodes, new_node]; 
        end

        % remove parent node
        nodes(1) = [];

        % Dont do too many iterations
        iter = iter +1;
        if iter>iterMax
            flag = -1;
        end
    end % while loop

    % Use control from the beginning of the history
    uInd = nodes(1).history(1);

end

% Output
U = params.hyd.Force2Torque(states(2))*params.hyd.ptoForceOptions(:);
out.controlValue = U(uInd);
out.controlIndex = uInd;
end



function node = initlizeNodes(nU,x0,params,ctrl,k,uPrev)
    node(nU) = struct();
    for uInd = 1:nU
        node(uInd).history(1) = uInd;
        [E_step,E_terminal,xf] = ComputeCost(params,ctrl,uInd,uPrev,x0,k,1);
        node(uInd).E_absorbed = E_step;
        node(uInd).cost = E_step + E_terminal;
        node(uInd).xf = xf;
    end
end


function [E_absorbed,E_terminal,xf] = ComputeCost(params,ctrl,uInd,uPrev,x,k,d)
% Comput the energy if we are picking a controller for block k and are
% anayling block d in the Astar algorithm.

U = params.hyd.Force2Torque(x(2))*params.hyd.ptoForceOptions(:);
u = U(uInd);

% Compute cost
% Energy in this block
w = ctrl.w;
b_exc = ctrl.b_exc(k+d-1);
Q_local = ctrl.Q_local;
E_mech = (w'*x+b_exc)*u + u^2*Q_local;
E_mech = (w'*x+b_exc)*u + u*ctrl.Q_local2*u/2;


% Switching Loss
E_sw = 0;

% Absorbed energy
E_absorbed = E_mech - E_sw;

% Advance state
M_block = ctrl.M_block;
T_block = ctrl.T_block(:,k+d-1);
H_block = ctrl.H_block;
xf = M_block*x+T_block+H_block*u;

% Future Cost
P = ctrl.termCost(d).P;
rT =ctrl.termCost(d).rT(k+d,:);
b = ctrl.termCost(d).b(k+d);
E_terminal = xf'*P*xf + rT*xf + b;

end