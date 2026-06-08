% This code is just a way of running main_WEC_Simulation.m with different
% parameters. The results from each set of parameters are collected here

clear, close all
%% Define cases to run
caseList = {
    'EHA',         'PI',             0;
    'EHA',         'MPC_QP',         0;
    'PassivePump', 'CoulombDamping', 2;
    'DHD',         'MPC_Astar_cont', 2;
    'DHD',         'MPC_Astar',      2;
    'DHD',         'MPC_Astar',      3;
    'DHD',         'MPC_Astar',      4;
    'DHD',         'PI',             2;
    'DHD',         'PI',             3;
    'DHD',         'PI',             4;
    'DHD',         'MPC_QP',         2;
    'DHD',         'MPC_QP',         3;
    'DHD',         'MPC_QP',         4;
    'DHD',         'Sliding Mode',   2;
    'DHD',         'Sliding Mode',   3;
    'DHD',         'Sliding Mode',   4;
};

simMatrix = cell2table(caseList, 'VariableNames', {'Drivetrain', 'Controller', 'PressureRails'});

%  Preallocate scalar columns for results
simMatrix.MechRGP = NaN(height(simMatrix), 1);
simMatrix.ElecRGP = NaN(height(simMatrix), 1);

% Preallocate cell columns to hold results that are matrices
% simMatrix.MatrixData = cell(height(simMatrix), 1);

%% 2. Loop Through and Populate Table
for i = 4%1:height(simMatrix)

    % make scturctur of parameters to pass around
    runParams = struct();
    runParams.drive = caseList{i, 1};
    runParams.controller = caseList{i, 2};
    runParams.pressure_rails = caseList{i, 3};

    % Display the current case being run
    disp('---------------------------------------------------------')
    if runParams.pressure_rails == 0
        fprintf('Running: %s with %s ...\n', runParams.drive, runParams.controller);
    else
        fprintf('Running: %d rail %s with %s ...\n', runParams.pressure_rails, runParams.drive, runParams.controller);
    end

    % Optimize the high pressure for the pressure rail cases
        % the medium rails are set to be evenly spaced
    switch runParams.drive
        case {'DHD','PassivePump'}
            % For optimizing pressure;
            % runParams.highPressure = OptPressure(runParams)*1e6;

            % For using a predetemrined value:
            runParams.highPressure = 35*1e6;
    end

    % Run the actual simulation
    main_WEC_Simulation;
    
    % Collect the results
    simMatrix.MechRGP(i) = eval.mechRGP;
    simMatrix.ElecRGP(i) = eval.elecRGP;
    
    % Assign full matrix into the cell column using curly braces {}
    % simMatrix.RawMatrixData{i} = eval.fullMatrix; 
end

%% Display Results
disp('--- Final Test Matrix ---');
disp(simMatrix); % Displays text and scalars cleanly

% How to extract the matrix from Row 4:
% row4Matrix = simMatrix.RawMatrixData{4};

% Generate latex code for this table:
% latexText = table2latex(simMatrix);

%% Other functions
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

function latexStr = table2latex(T, varargin)
% TABLE2LATEX Converts a MATLAB table into a LaTeX tabular text string.

    % Set up optional precision format argument
    p = inputParser;
    addParameter(p, 'Precision', "%.4f", @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    floatFmt = string(p.Results.Precision);

    % Extract structural info
    varNames = string(T.Properties.VariableNames);
    [numRows, numCols] = size(T);
    
    % Initialize string builder array (Using modern string array)
    lines = string.empty;
    
    % 1. Setup table alignment
    colAlignment = repmat("c", 1, numCols).join(""); 
    lines(end+1) = sprintf("\\begin{tabular}{%s}", colAlignment);
    lines(end+1) = "\hline";
    
    % 2. Generate column headers
    headerRow = "";
    for j = 1:numCols
        % Escape underscores in headers if present
        cleanHeader = strrep(varNames(j), "_", "\_");
        headerRow = headerRow + sprintf("\\textbf{%s}", cleanHeader);
        if j < numCols
            headerRow = headerRow + " & ";
        end
    end
    lines(end+1) = headerRow + " \\\\";
    lines(end+1) = "\hline";
    
    % 3. Loop through table rows and values
    for i = 1:numRows
        rowStr = "";
        for j = 1:numCols
            val = T{i, j};
            
            % Handle nested cells
            if iscell(val)
                val = val{1};
            end
            
            % Format based on data type
            if isnumeric(val)
                if isnan(val)
                    cellStr = "-";
                elseif mod(val, 1) == 0
                    cellStr = sprintf("%d", val);
                else
                    cellStr = sprintf(floatFmt, val);
                end
            else
                % Convert text and safely escape LaTeX underscores (e.g., MPC_QP -> MPC\_QP)
                cellStr = strrep(string(val), "_", "\_");
            end
            
            rowStr = rowStr + cellStr;
            if j < numCols
                rowStr = rowStr + " & ";
            end
        end
        lines(end+1) = rowStr + " \\\\";
    end
    
    % 4. Close table
    lines(end+1) = "\hline";
    lines(end+1) = "\end{tabular}";
    
    % Combine lines using a clean newline character
    latexStr = join(lines, newline);
    
    % Print to command window for quick copy-pasting
    fprintf('\n--- Copy the text below ---\n\n%s\n\n---------------------------\n', latexStr);
end

