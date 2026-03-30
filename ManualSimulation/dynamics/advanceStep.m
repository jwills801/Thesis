function finalState = advanceStep(initialState,dt,sys,T)

% If the position is near the end stops, apply an end stop force
pos = initialState(2);
vel = initialState(1);
if abs(pos) > pi/4
    endStopTorque = -pos*1e7 - vel*1e7;
else
    endStopTorque = 0;
end

% Calculate xdot
xDot = sys.A*initialState + sys.B*(T+endStopTorque) ;

% Forward euler one step
    finalState = initialState + xDot*dt;
    
end