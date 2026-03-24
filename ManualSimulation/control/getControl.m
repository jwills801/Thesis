function cntrl = getControl(params,wave)

% controller = 'PI';
controller = 'Sliding Mode';
    cntrl.timeHorizon = .2; % Amount of time to simulate
    cntrl.horizonInd = round(cntrl.timeHorizon/params.simu.dt);
    cntrl.lambda = 1; % Defines the sliding surface
    cntrl.phi = 3e-2; % band around sliding surface

optTraj = getOptimal(params,wave);

cntrl.optTraj = optTraj;
cntrl.controller = controller;

end