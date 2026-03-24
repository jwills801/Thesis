function [u, uInd] = controlLaw(params,cntrl,wave,states,uInd_history)

switch cntrl.controller
    case 'PI'
        out = PIcontrol(params,states,uInd_history);
    case 'Sliding Mode'
        out = slidingMode(params,cntrl,wave,states,uInd_history);
end
u = out.controlValue;
uInd = out.controlIndex;
end