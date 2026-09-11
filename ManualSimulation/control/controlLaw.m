% controlLaw.m
% Dispatches to the controller named in params.runParams.controller and
% returns its chosen control torque for the current timestep.
% Calls: PIcontrol.m, slidingMode.m, MPC_QP.m, MPC_Astar.m,
%   MPC_Astar_cont.m, coulombDamping.m
% Called by: dynamics/timeLoop.m, diagnostics/checkAstarVsBruteForce.m
function [u, uInd] = controlLaw(params,ctrl,wave,states,uInd_history)

switch params.runParams.controller
    case 'PI'
        out = PIcontrol(params,ctrl,states,uInd_history);
    case 'SlidingMode'
        out = slidingMode(params,ctrl,wave,states,uInd_history);
    case 'MPC_QP'
        out = MPC_QP(params,ctrl,wave,states,uInd_history);
    case 'MPC_Astar'
        out = MPC_Astar(params,ctrl,wave,states,uInd_history);
    case 'MPC_Astar_cont'
        out = MPC_Astar_cont(params,ctrl,wave,states,uInd_history);
    case 'CoulombDamping'
        out = coulombDamping(params,states);
end
u = out.controlValue;
uInd = out.controlIndex;
end