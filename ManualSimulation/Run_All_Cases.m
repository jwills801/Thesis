clear, close all
%% 1. Define the Master Test Matrix
caseList = {
    'PassivePump', 'CoulombDamping', 0;
    'EHA',         'PI',             0;
    'EHA',         'MPC_QP',         0;
    'DHD',         'PI',             2;
    'DHD',         'PI',             3;
    'DHD',         'PI',             4;
    'DHD',         'MPC_Astar',      2;
    'DHD',         'MPC_Astar',      3;
    'DHD',         'MPC_Astar',      4;
    'DHD',         'MPC_QP',         2;
    'DHD',         'MPC_QP',         3;
    'DHD',         'MPC_QP',         4;
};

% Create the base table
simMatrix = cell2table(caseList, 'VariableNames', {'Drivetrain', 'Controller', 'PressureRails'});

%  Preallocate scalar columns
simMatrix.MechRGP = NaN(height(simMatrix), 1);
simMatrix.ElecRGP = NaN(height(simMatrix), 1);

% Preallocate cell columns to hold full matrices of varying sizes
% simMatrix.RawMatrixData = cell(height(simMatrix), 1);

%% 2. Loop Through and Populate Table
for i = 8%1:height(simMatrix)
    runParams = struct();
    runParams.drive = caseList{i, 1};
    runParams.controller = caseList{i, 2};
    % Extract inputs
    drivetrain     = simMatrix.Drivetrain{i};
    controller     = simMatrix.Controller{i};
    pressure_rails = simMatrix.PressureRails(i);
    
        if pressure_rails == 2
            runParams.pressureRails = [0 35]*1e6;
        elseif pressure_rails == 3
            runParams.pressureRails = [0 17 35]*1e6;
        elseif pressure_rails == 4
            runParams.pressureRails = [0 8 27 35]*1e6;
        elseif pressure_rails == 5
            runParams.pressureRails = [0 3.5 18.5 29 35]*1e6;
        end

        disp('---------------------------------------------------------')
        if pressure_rails == 0
            fprintf('Running: %s with %s ...\n', runParams.drive, runParams.controller);
        else
            fprintf('Running: %d rail %s with %s ...\n', pressure_rails, runParams.drive, runParams.controller);
        end
        
    main_WEC_Simulation;
    
    % Assign scalars
    simMatrix.MechRGP(i) = eval.mechRGP;
    simMatrix.ElecRGP(i) = eval.elecRGP;
    
    % Assign full matrix into the cell column using curly braces {}
    % simMatrix.RawMatrixData{i} = eval.fullMatrix; 
end

%% 3. Display and Access Data
disp('--- Final Test Matrix ---');
disp(simMatrix); % Displays text and scalars cleanly

% How to extract the matrix from Row 4:
% row4Matrix = simMatrix.RawMatrixData{4};

