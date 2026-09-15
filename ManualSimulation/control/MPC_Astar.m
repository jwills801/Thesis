% MPC_Astar.m
% DHD's receding-horizon A* branch-and-bound search over discrete PTO
% force combinations (cap-rail x rod-rail), using ctrl.termCost (from
% getControl.m) as the terminal-cost heuristic and, if
% params.runParams.considerLosses, the switch-loss map for switching
% costs. Self-contained (local initlizeNodes/ComputeCost/getSwitchingLoss).
% out.capHit is true when this window's search was truncated by
% astarIterMax before reaching its normal d==m_Astar exit (also raises a
% MPC_Astar:iterCapHit warning); controlLaw.m/timeLoop.m/evaluate.m
% propagate this through to eval.nAstarCapHits for the whole simulation.
% Calls: none
% Called by: control/controlLaw.m
function out = MPC_Astar(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control input indexes
out.capHit = false; % overwritten below only if this window's search hits astarIterMax
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
    % iterMax overridable via runParams.astarIterMax (default 1000,
    % unchanged for every existing caller). The node list grows by
    % (nU-1) per iteration and is fully re-sorted every iteration (an
    % O(n^2) cost, not fixed here -- see diagnostics/ReadMe.md), which is
    % mild for small nU (e.g. 2-rail DHD, nU=4) but severe for large nU
    % (4-rail DHD, nU=16: node list reaches ~15000 entries at iter=1000,
    % re-sorted every one of ~2500 control windows per simulation --
    % observed to take many hours for a single sim). Capping iterMax
    % lower bounds worst-case runtime at the cost of the search
    % sometimes falling back before reaching its normal d==m_Astar exit,
    % i.e. a smaller compute budget, not a change to the search itself.
    flag = 0; iter = 0;
    if isfield(params.runParams,'astarIterMax')
        iterMax = params.runParams.astarIterMax;
    else
        iterMax = 1e3;
    end
    while flag == 0

        % sort nodes
        [~, idx] = sort([nodes.cost], 'ascend');

        % reorder nodes
        nodes = nodes(idx);

        % Get info about this node
        d = length(nodes(1).history);
        uIndPrev = nodes(1).history(end);
        xf = nodes(1).xf;

        % Check if we reached the end of the horizon
        if d == ctrl.m_Astar
            flag = 1;
            break
        end

        % look at all the possible options from this node
        for uInd = 1:nU
            [E_step, E_term, x_next] = ComputeCost(params, ctrl, uInd, uIndPrev, xf, k, d+1);
        
            new_node.history = [nodes(1).history, uInd];
            new_node.xf = x_next;
            new_node.E_absorbed = nodes(1).E_absorbed + E_step;
            new_node.cost = new_node.E_absorbed + E_term;

            % Add child to the pool
            nodes = [nodes, new_node]; 
        end

        % remove parent node
        nodes(1) = [];

        % Dont do too many iterations
        iter = iter +1;
        if iter>iterMax
            flag = -1;
            out.capHit = true;
            % A cap this low should be rare/never for a well-configured
            % run -- see astarIterMax's callers for the intended budget.
            % If this fires often within one simulation, the search is
            % being truncated before it naturally converges (d==m_Astar),
            % which silently degrades this window's chosen control -- not
            % a crash, so it needs this warning to be noticed at all.
            % timeLoop.m/evaluate.m aggregate these into
            % eval.nAstarCapHits so a whole simulation's total is visible
            % without grepping logs.
            warning('MPC_Astar:iterCapHit', ...
                ['A* search hit its %d-iteration cap at t=%.3fs (m_Astar=%d, ' ...
                 'nU=%d) before reaching its normal d==m_Astar exit.'], ...
                iterMax, params.simu.time(fineTimeInd), ctrl.m_Astar, nU);
        end
    end % while loop

    % Use control from the beginning of the history
    uInd = nodes(1).history(1);

    % Efficiency
    fullTree = (nU^ctrl.m_Astar-1)/(nU-1);
    AStar_eff = iter/fullTree; %#ok<NASGU>

    if ~isfinite(nodes(1).cost)
        warning('MPC_Astar:nonfiniteCost', ...
            'A* chose a non-finite-cost node at t=%.3fs.', params.simu.time(fineTimeInd));
    end

end

% Output
U = params.hyd.Force2Torque(states(2))*params.hyd.ptoForceOptions(:);
out.controlValue = U(uInd);
out.controlIndex = uInd;
end



function node = initlizeNodes(nU,x0,params,ctrl,k,uIndPrev)
    node(nU) = struct();
    for uInd = 1:nU
        node(uInd).history(1) = uInd;
        [E_step,E_terminal,xf] = ComputeCost(params,ctrl,uInd,uIndPrev,x0,k,1);
        node(uInd).E_absorbed = E_step;
        node(uInd).cost = E_step + E_terminal;
        node(uInd).xf = xf;
    end
end


function [E_absorbed,E_terminal,xf] = ComputeCost(params,ctrl,uInd,uIndPrev,x,k,d)
% Compute the energy if we are picking a controller for block k and are
% anayling block d in the Astar algorithm.

u = params.hyd.Force2Torque(x(2))*params.hyd.ptoForceOptions(uInd);

% Compute cost
% Energy in this block
w = ctrl.w;
b_exc = ctrl.b_exc(k+d-1);
Q_local = ctrl.Q_local;
E_mech = (w'*x+b_exc)*u + u^2*Q_local;

% Switching Loss
if params.runParams.considerLosses
    E_sw = getSwitchingLoss(params,x,uInd,uIndPrev);
else
    E_sw = 0;
end

% Absorbed energy
E_absorbed = E_mech + E_sw;

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


function E_sw = getSwitchingLoss(params,x,uInd,uIndPrev)
switchMap = params.hyd.switchMap;

% Find cap and rod side volumes and flows
[cap,rod] = params.hyd.getVolandFlow(params,x);


% Which pressure rail did we switch from?
    % params.hyd.ptoTorqueOptions is a matrix
    % Each row is a different cap side option
    % Each col is a different rod side option
[cap.switchFromInd,rod.switchFromInd] = ind2sub(size(params.hyd.ptoForceOptions),uIndPrev);
[cap.switchToInd,rod.switchToInd] = ind2sub(size(params.hyd.ptoForceOptions),uInd);

% compute loss on each side
E_sw = 0;
for side = [cap, rod]
    % define variable to be interpolated on
    switchFrom = params.hyd.pressureRails(side.switchFromInd);
    switchTo = params.hyd.pressureRails(side.switchToInd);
    switchVelA = side.velA;
    switchVol = side.vol + switchMap.hoseVolume;

    % interpolate. With a large (or no) astarIterMax, the search's own
    % multi-step lookahead (ComputeCost propagating xf several steps
    % under a fixed candidate force) can predict velA/vol well outside
    % switchMap's fixed +-1.5*capArea/[hoseVolume,hoseVolume+stroke*
    % capArea] range -- a hypothetical branch several steps deep under a
    % sustained aggressive rail choice has no reason to stay in-range,
    % unlike the true simulated trajectory. Passing an explicit
    % extrapolation value makes any such branch cost a fixed large
    % penalty instead of silently NaN (observed: 89 non-finite-cost
    % windows in one DHD2 run once astarIterMax stopped truncating the
    % search before it reached this regime) -- finite so it still sorts
    % correctly (unlike NaN) and sums cleanly with the rest of the cost,
    % but large enough that a branch this far outside the map's modeled
    % range never wins a comparison against an in-range one.
    OUT_OF_RANGE_PENALTY = 1e12; % J -- swamps any real switching loss (O(1e2)-O(1e6) J)
    loss = interpn(switchMap.PR,switchMap.PR,switchMap.velA_vals,switchMap.vol_vals, switchMap.Eloss,...
        switchFrom,switchTo,switchVelA,switchVol,'linear',OUT_OF_RANGE_PENALTY);

    % Add up loss from each side
    E_sw = E_sw + loss;
end

end