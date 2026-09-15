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
% See README.md (repo root) for the full folder structure and data-flow
% diagram -- not duplicated here to avoid the two going stale independently.

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
addpath("models/")
addpath("optimization/")
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