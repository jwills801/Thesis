% timeLoop.m
% The main time-domain simulation loop: at every fine timestep, asks
% controlLaw.m for a control torque and advances the true plant one step
% (advanceStep.m). Returns the full state/control/time history in dyn.
% Calls: control/controlLaw.m, dynamics/advanceStep.m
% Called by: main_WEC_Simulation.m, parameters/optimizePressure.m,
%   parameters/sizeCylinderArea.m, diagnostics/checkEnergyBalance.m,
%   validatePhase2Subset.m
% (checkAstarVsBruteForce.m calls controlLaw.m directly instead, to
% inspect a single decision without running the full time loop.)
function dyn = timeLoop(params,wave,ctrl)
% Time
t = params.simu.time;
dt = params.simu.dt;
sys = params.phys.sys;
Texc = wave.torque.Texc;
    
% initilize state vector and set I.C.
states = NaN(length(sys.A),length(t));
states(:,1) = zeros(length(sys.A),1);

% Initilize control and set I.C.
uInd = ones(length(t),1); uInd(1) = 1; 
u = NaN(length(t),1); u(1) = 0;

% This is the main time-domain simulation (can be tens of thousands of
% steps), so it's worth a waitbar -- but only update it every ~1% of
% progress rather than every single step.
waitbarObj = waitbar(0,'Simulating WEC Dynamics');
updateEvery = max(1,round(length(t)/100));
for timeInd = 1:length(t)-1
    if mod(timeInd,updateEvery) == 0
        waitbar(timeInd/length(t),waitbarObj);
    end

    [u(timeInd), uInd(timeInd)] = controlLaw(params,ctrl,wave,states(:,timeInd),uInd(1:timeInd-1));

    states(:,timeInd+1) = advanceStep(states(:,timeInd),dt,sys,u(timeInd)+Texc(timeInd),params.phys.b);
end
close(waitbarObj)

% For the last time step
u(timeInd+1) = u(timeInd);

% output results
dyn.u = u;
dyn.uInd = uInd;
dyn.states = states;
dyn.thetaDot = states(1,:)';
dyn.theta = states(2,:)';
dyn.t = t;

end