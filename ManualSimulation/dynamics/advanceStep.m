function finalState = advanceStep(initialState,dt,sys,T)

% Forward euler one step
    finalState = initialState + (sys.A*initialState + sys.B*T)*dt;
    
% lsim (takes many steps)
    % statesFine = lsim(sys,[T T],[0 dt],initialState);
    % finalState = statesFine(end,:);
end