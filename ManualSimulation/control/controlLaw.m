function [u, uInd] = controlLaw(params,ctrl,wave,states,uInd_history)

switch ctrl.controller
    case 'PI'
        out = PIcontrol(params,ctrl,states,uInd_history);
    case 'Sliding Mode'
        out = slidingMode(params,ctrl,wave,states,uInd_history);
            case 'MPC'
        out = MPC(params,ctrl,wave,states,uInd_history);
    case 'Coulomb Damping'
        out = coulombDamping(params,states);
end
u = out.controlValue;
uInd = out.controlIndex;
end