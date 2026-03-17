function dyn = timeLoop(params,wave,optTraj)
% Time
t = params.simu.time;
dt = params.simu.dt;
sys = params.phys.sys;
Texc = wave.torque.Texc;
    
% initilize state vector and set I.C.
states = NaN(length(sys.A),length(t));
states(:,1) = zeros(length(sys.A),1);

% Initilize control and set I.C.
u = zeros(1,length(t)); u(1) = 0;

waitbarObj = waitbar(0,'Simulating WEC Dynamics');
for timeInd = 1:length(t)-1
    waitbar(timeInd/length(t),waitbarObj);

    u(timeInd) = PIcontrol(states(:,timeInd),params);

    states(:,timeInd+1) = advanceStep(states(:,timeInd),dt,sys,u(timeInd)+Texc(timeInd));
end
close(waitbarObj)

% For the last time step
u(timeInd+1) = u(timeInd);

% output results
dyn.u = u;
dyn.states = states;
dyn.thetaDot = states(2,:);
dyn.theta = states(2,:);
dyn.t = t;

end