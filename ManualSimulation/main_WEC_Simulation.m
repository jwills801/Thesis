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
        % getOptimal.m
        % PIcontrol.m
        % slidingMode.m
    % dynamics/
        % timeLoop.m
        % advanceStep.m
    % evaluation/                           # Package for loss analysis
        % analyzeValves.m
        % calculateKinematics.m
        % findEvents.m
    % plotting/                       # Package for visualization
        % plotAll.m
        % plotControl.m
        % plotLosses.m


clear, close all

% Load parameters
addpath("parameters/")
params = getParameters;

% Calculate excitation torque
addpath("wave/")
wave = generateExcitingTorque(params);

% Initialize control
addpath("control/")
optTraj = getOptimal(params,wave);

% Simulate Dynamics
addpath("dynamics/")
dyn = timeLoop(params,wave,optTraj);

% Evaluate
addpath("evaluation/")
eval = evaluate(params,dyn);

%% Plot
addpath("plotting/")
plotting(params,wave,optTraj,dyn,eval);

