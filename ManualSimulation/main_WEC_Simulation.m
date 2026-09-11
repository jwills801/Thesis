% main_WEC_Simulation.m
% Top-level script: given runParams (drive/controller/pressure_rails/
% considerLosses etc. already in the workspace), sets drivetrain-specific
% defaults, builds params/wave/ctrl, runs the time loop, evaluates, and
% plots. Run directly (sets its own runParams first) or via Run_All_Cases.m
% (which sets runParams per case then calls this as a script).
% Calls: parameters/getParameters.m (-> getHydraulic.m, getPhysical.m,
%   getSimulation.m), wave/generateExcitingTorque.m, control/getControl.m,
%   dynamics/timeLoop.m, evaluation/evaluate.m, plotting/plotAll.m
% Called by: Run_All_Cases.m (as a script, once per case row)
%
%% Code Structure
% main_WEC_Simulation.m           # Top-level script
    % parameters/                     # Package for parameter functions
        % getHydraulic.m
        % getParameters.m
        % getPhysical.m
        % getSimulation.m
        % makeEHALossMap.m
        % makeSwitchLossMap.m
    % wave/                           # Package for wave functions
        % calculateWavePower.m
        % generateExcitingTorque.m
        % getSpectrum.m	
        % getTorqueTimeSeries.m
    %control/                        # Package for control functions
        % controlLaw.m
        % coulombDamping.m
        % getControl.m
        % getOptimal.m
        % MPC_Astar_cont.m
        % MPC_Astar.m
        % MPC_QP.m
        % PIcontrol.m
        % slidingMode.m
    % dynamics/
        % timeLoop.m
        % advanceStep.m
    % evaluation/                      # Package for loss analysis
        % evaluate.m
        % getEHALoss.m
        % getHECMLoss.m
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
        % Respect a pre-set considerLosses (mechanical- vs electrical-
        % energy-optimized EHA, see Run_All_Cases.m's ConsiderLosses
        % column); default to mechanical-optimized if not specified.
        if ~isfield(runParams,'considerLosses')
            runParams.considerLosses = 0;
        end
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
%%

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
eval = evaluate(params,dyn,ctrl);

% Plot
addpath("plotting/")
eval = plotAll(params,wave,ctrl,dyn,eval);

toc