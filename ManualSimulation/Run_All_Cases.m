clear, close all
% Define the master list of cases
caseList = {
    'PassivePump', 'CoulombDamping';
    'EHA',         'PI';
    'EHA',         'MPC_QP';
    'DHD',         'PI';
    'DHD',         'MPC_Astar';
    'DHD',         'MPC_QP'
};

% Initialize results structure if it doesn't exist
eval_results = struct();

for i = 4%:size(caseList, 1)
    runParams = struct();
    runParams.drive = caseList{i, 1};
    runParams.controller = caseList{i, 2};
    
    % Determine pressure rails to run
    if strcmp(runParams.drive, 'DHD')
        railsToRun = [2, 3, 4];
        railsToRun = 3;
    elseif strcmp(runParams.drive, 'PassivePump')
        railsToRun = 2;
    else
        railsToRun = 0; % Dummy value for non-DHD cases
    end
    
    for p = railsToRun
        % Update simulation parameters
        % (Assuming your model uses these workspace variables)
        if p == 2
            runParams.pressureRails = [0 35]*1e6;
        elseif p == 3
            runParams.pressureRails = [0 17 35]*1e6;
        elseif p == 4
            runParams.pressureRails = [0 17 22 35]*1e6;
        elseif p == 5
            runParams.pressureRails = [0 8 17 26 35]*1e6;
        end
        disp('---------------------------------------------------------')
        fprintf('Running: %s with %s (Rails: %d)...\n', runParams.drive, runParams.controller, p);
        
        main_WEC_Simulation;

        % Generate a valid field name for the structure
        % e.g., DHD_MPC_QP_3rails
        caseName = sprintf('%s_%s', runParams.drive, runParams.controller);
        if p > 0
            caseName = sprintf('%s_%drails', caseName, p);
        end
        
        % Store results in the eval structure
        eval_results.(caseName) = eval; 
        eval_results.(caseName).status = 'Completed'; % Placeholder
    end
end

disp('All simulation cases processed.');

%% Make a table from the results
% Get all case names from the structure
caseNames = fieldnames(eval_results);

% Preallocate cell arrays for the table columns
CaseName = caseNames;
MechRGP = zeros(length(caseNames), 1);
ElecRGP = zeros(length(caseNames), 1);

% Extract values from each case
for i = 1:length(caseNames)
    currentCase = caseNames{i};
    
    % Access the nested eval structure for the current case
    % Note: Adjusting the path to eval_results.(currentCase).eval 
    % based on your typical simulation data logging
    MechRGP(i) = eval_results.(currentCase).mechRGP;
    ElecRGP(i) = eval_results.(currentCase).elecRGP;
end

% Create the final table
resultsTable = table(CaseName, MechRGP, ElecRGP);

% Display the table in the Command Window
disp(resultsTable);

