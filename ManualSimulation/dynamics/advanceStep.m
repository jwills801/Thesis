% advanceStep.m
% Advances the flap's true plant state one forward-Euler step under the
% linear state-space system sys, plus an end-stop torque (large stiffness
% if |theta|>pi/4) and optional quadratic viscous damping b|thetaDot|thetaDot.
% This is the TRUE plant; controllers' internal prediction models
% (getTransition.m) use ZOH instead, and stay on the zero-damping linear
% system regardless of b.
% Calls: none
% Called by: dynamics/timeLoop.m
function finalState = advanceStep(initialState,dt,sys,T,b)
% b: quadratic viscous damping coefficient (main(1).tex tau_v=b|thetaDot|thetaDot),
% applied here to the true plant only. Optional, defaults to 0 (no damping,
% matches every existing result) so existing call sites are unaffected.
if nargin < 5
    b = 0;
end

% If the position is near the end stops, apply an end stop force
pos = initialState(2);
vel = initialState(1);
if abs(pos) > pi/4
    endStopTorque = -pos*1e7 - vel*1e7;
else
    endStopTorque = 0;
end

% Quadratic viscous damping torque, opposing motion
viscousTorque = -b*abs(vel)*vel;

% Calculate xdot
xDot = sys.A*initialState + sys.B*(T+endStopTorque+viscousTorque) ;

% Forward euler one step
    finalState = initialState + xDot*dt;
    
end