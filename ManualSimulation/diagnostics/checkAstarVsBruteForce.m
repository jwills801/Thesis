% checkAstarVsBruteForce.m
% Confirms MPC_Astar.m's branch-and-bound search returns the TRUE minimum
% cost sequence, by brute-force enumerating every control sequence over a
% short 2-step horizon (ctrl.m_Astar reduced to 2) and comparing the
% chosen first-step control index. This validates the search logic itself,
% not just that the cost formulas are right (checked separately).
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, controlLaw.m (which dispatches to MPC_Astar.m;
%   not called directly -- its own stepCost/switchingLoss sub-functions
%   duplicate MPC_Astar.m's cost formulas independently for comparison)
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'wave'));
addpath(fullfile(root,'control'));

runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',2, ...
    'considerLosses',1,'rodArea',(0.0254*6)^2*pi,'capArea',1.5*(0.0254*6)^2*pi, ...
    'highPressure',35e6);
params = getParameters(runParams);
wave = generateExcitingTorque(params);
ctrl = getControl(params,wave); % full production-scale precompute
ctrl.m_Astar = 2;                % shrink search horizon for exhaustive brute force

x0 = zeros(4,1); % matches timeLoop.m's initial condition
uIndPrev = 1; k = 1;
nU = length(params.hyd.ptoForceOptions(:));

[~,uInd_astar] = controlLaw(params,ctrl,wave,x0,[]);

bestCost = inf; bestFirst = NaN;
for u1 = 1:nU
    [E1,xf1] = stepCost(params,ctrl,u1,uIndPrev,x0,k,1);
    for u2 = 1:nU
        [E2,xf2] = stepCost(params,ctrl,u2,u1,xf1,k,2);
        P = ctrl.termCost(2).P; rT = ctrl.termCost(2).rT(k+2,:); b = ctrl.termCost(2).b(k+2);
        cost = E1+E2 + xf2'*P*xf2 + rT*xf2 + b;
        if cost < bestCost
            bestCost = cost; bestFirst = u1;
        end
    end
end

fprintf('--- checkAstarVsBruteForce ---\n');
if uInd_astar == bestFirst
    fprintf('PASS  A* first-step choice (%d) matches brute-force optimum (%d)\n',uInd_astar,bestFirst);
else
    fprintf('FAIL  A* first-step choice (%d) does NOT match brute-force optimum (%d)\n',uInd_astar,bestFirst);
end

function [E,xf] = stepCost(params,ctrl,uInd,uIndPrev,x,k,d)
u = params.hyd.Force2Torque(x(2))*params.hyd.ptoForceOptions(uInd);
w = ctrl.w; b_exc = ctrl.b_exc(k+d-1); Q_local = ctrl.Q_local;
E_mech = (w'*x+b_exc)*u + u^2*Q_local;
if params.runParams.considerLosses
    E_sw = switchingLoss(params,x,uInd,uIndPrev);
else
    E_sw = 0;
end
E = E_mech + E_sw;
xf = ctrl.M_block*x + ctrl.T_block(:,k+d-1) + ctrl.H_block*u;
end

function E_sw = switchingLoss(params,x,uInd,uIndPrev)
switchMap = params.hyd.switchMap;
[cap,rod] = params.hyd.getVolandFlow(params,x);
[cap.switchFromInd,rod.switchFromInd] = ind2sub(size(params.hyd.ptoForceOptions),uIndPrev);
[cap.switchToInd,rod.switchToInd] = ind2sub(size(params.hyd.ptoForceOptions),uInd);
E_sw = 0;
for side = [cap, rod]
    switchFrom = params.hyd.pressureRails(side.switchFromInd);
    switchTo = params.hyd.pressureRails(side.switchToInd);
    switchVelA = side.velA; switchVol = side.vol + switchMap.hoseVolume;
    loss = interpn(switchMap.PR,switchMap.PR,switchMap.velA_vals,switchMap.vol_vals,switchMap.Eloss,...
        switchFrom,switchTo,switchVelA,switchVol);
    E_sw = E_sw + loss;
end
end
