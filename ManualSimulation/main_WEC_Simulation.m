%% Code Structure
% main_WEC_Simulation.m           # Top-level script
    % parameters/                     # Package for parameter functions
        % getParameters.m
        % getPhysical.m
        % getHydraulic.m
        % getSimulation.m
        % getControl.m
    % wave/                           # Package for wave functions
        % generateExcitingTorque.m
        % getSpectrum.m	
        % getTorqueTimeSeries.m
        % calculateWavePower.m
    %control/                        # Package for control functions
        % getControl.m
        % controlLaw.m
        % getOptimal.m
        % coulombDamping.m
        % PIcontrol.m
        % slidingMode.m
    	% MPC_DP.m
        % MPC.m
    % dynamics/
        % timeLoop.m
        % advanceStep.m
    % evaluation/                           # Package for loss analysis
        % evaluate.m
        % makeSwitchLossMap.m
        % getValveLoss.m
    % plotting/                       # Package for visualization
        % plotAll.m

% Assign parameters based on drivetrain selection
switch runParams.drive
    case 'PassivePump'
        runParams.rodArea = (0.0254*6)^2*pi;
        runParams.capArea = 1.5*runParams.rodArea;
        runParams.highPressure = 20.6*1e6;
    case 'EHA'
        runParams.considerLosses = 0;
        runParams.rodArea = (0.0254*8)^2*pi;
        runParams.capArea = runParams.rodArea;
    case 'DHD'
        runParams.considerLosses = 1;
        runParams.rodArea = (.0254*6)^2*pi; % m^2: Radius squared times pi
        runParams.capArea = 1.5*runParams.rodArea; % m^2: Area ratio times rod Area
        runParams.highPressure = 33.3*1e6;
end



%%
tic
% Load parameters
addpath("parameters/")
params = getParameters(runParams);

% Calculate excitation torque
addpath("wave/")
wave = generateExcitingTorque(params);

% Initialize control
addpath("control/")
ctrl = getControl(params,wave);

%% Simulate Dynamics
addpath("dynamics/")
dyn = timeLoop(params,wave,ctrl);

%% Evaluate
addpath("evaluation/")
eval = evaluate(params,dyn);

% Plot
addpath("plotting/")
eval = plotAll(params,wave,ctrl,dyn,eval);

toc