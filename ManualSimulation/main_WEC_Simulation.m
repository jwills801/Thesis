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

clear, close all

%% Run defining Parameters
runParams = struct();

% Select Drivetrain and Controller
% drivetrain = 'PassivePump'; controller = 'CoulombDamping';
drivetrain = 'EHA'; controller = 'PI';
% drivetrain = 'EHA'; controller = 'MPC_QP';
% drivetrain = 'DHD'; controller = 'PI';
% drivetrain = 'DHD'; controller = 'MPC_Astar';
runParams.drive = drivetrain;
runParams.controller = controller;

% Assign parameters based on drivetrain selection
switch drivetrain
    case 'PassivePump'
        runParams.pressureRails = [0 35e6];
        runParams.rodArea = (0.0254*2)^2*pi;
        runParams.capArea = 1.5*runParams.rodArea;
    case 'EHA'
        runParams.rodArea = (0.0254*8)^2*pi;
        runParams.capArea = runParams.rodArea;
    case 'DHD'
        runParams.considerSwitchingLoss = 1;
        runParams.pressureRails = [0 10 20 35]*1e6;
        runParams.rodArea = (.0254*6)^2*pi; % m^2: Radius squared times pi
        runParams.capArea = 1.5*runParams.rodArea; % m^2: Area ratio times rod Area
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

% Simulate Dynamics
addpath("dynamics/")
dyn = timeLoop(params,wave,ctrl);

% Evaluate
addpath("evaluation/")
eval = evaluate(params,dyn);

%% Plot
addpath("plotting/")
plotAll(params,wave,ctrl,dyn,eval);

toc