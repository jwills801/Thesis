function out = MPC_Astar_cont2(params,ctrl,wave,states,uInd_history)
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

    uInd = params.uInd(k);

end

% Output
out.controlValue = getTorque(params,states,uInd);
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

function u = getTorque(params,states,uInd)
% Find the current force rail
u_tmp = params.hyd.Force2Torque(states(2))*params.hyd.ptoForceOptions(uInd);

% now subtact the force that the electric motor would take out
u = u_tmp - params.damping*states(1);
end


function [E_absorbed,E_terminal,xf] = ComputeCost(params,ctrl,uInd,uIndPrev,x,k,d)
% Compute the energy if we are picking a controller for block k and are
% anayling block d in the Astar algorithm.

u = getTorque(params,x,uInd);

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

    % interpolate
    loss = interpn(switchMap.PR,switchMap.PR,switchMap.velA_vals,switchMap.vol_vals, switchMap.Eloss,...
        switchFrom,switchTo,switchVelA,switchVol);

    % Add up loss from each side
    E_sw = E_sw + loss;
end

end