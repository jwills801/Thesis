clear, close all
%% 1. Define the Master Test Matrix
caseList = {
    'PassivePump', 'CoulombDamping', 2;
    'EHA',         'PI',             0;
    'EHA',         'MPC_QP',         0;
    'DHD',         'PI',             2;
    'DHD',         'PI',             3;
    'DHD',         'PI',             4;
    'DHD',         'MPC_Astar_cont', 2;
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
for i = 7%1:height(simMatrix)
    runParams = struct();
    runParams.drive = caseList{i, 1};
    runParams.controller = caseList{i, 2};
    runParams.pressure_rails = caseList{i, 3};

        disp('---------------------------------------------------------')
        if runParams.pressure_rails == 0
            fprintf('Running: %s with %s ...\n', runParams.drive, runParams.controller);
        else
            fprintf('Running: %d rail %s with %s ...\n', runParams.pressure_rails, runParams.drive, runParams.controller);
        end
switch runParams.drive
    case {'DHD','PassivePump'}
        % runParams.highPressure = OptPressure(runParams)*1e6;
    %case 
    %    runParams.highPressure = 35*1e6;
end
% runParams.highPressure = 35*1e6;
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

function P  = OptPressure(runParams)
    iter = 0; itermax =10;
    flag = 0;
    P = 30; Pprev = 0; Jprev = 0;

    while flag ==0
    runParams.highPressure = P*1e6;
    main_WEC_Simulation; close all
    J = eval.elecRGP*100;

    Pnew = P + 5*(J-Jprev)/(P-Pprev)

    % If we have converged, then leave the loop
    if (abs(P-Pprev) < 1) && (abs(J-Jprev) < 1)
        flag = 1;
    end

    % update for next cycle
    Jprev = J; Pprev = P; P = Pnew;

    % If we have too many iterations then leave
    iter = iter +1;
    if iter > itermax
        flag = -1;
        disp('Max iterations reached on pressure optimization')
    end
    end
end

