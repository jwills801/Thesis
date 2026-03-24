function out = coulombDamping(params,states)
thetaDot = states(1);
theta = states(2);

ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);


if thetaDot >=0
    [u,uInd] = min(ptoTorqueOptions);
else
    [u,uInd] = max(ptoTorqueOptions);
end

% Smooth with hyberbolic tangent
u_pos = min(ptoTorqueOptions);
u_neg = max(ptoTorqueOptions);

smooth_sgn = tanh(thetaDot / .01);
u = ( (u_pos - u_neg)/2 + (u_pos + u_neg)/2*smooth_sgn ) *smooth_sgn;

% output results
out.controlValue = u;
out.controlIndex = uInd;