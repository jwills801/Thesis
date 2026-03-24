function out = slidingMode(params,ctrl,wave,states,uInd_history)
% uInd_history is a vectory of the previous control inputs

% Calculate the amount of time since the last switch
timeInd = length(uInd_history)+1;
lastSwitchInd = find(diff(uInd_history)~=0,1,'last');
    if isempty(lastSwitchInd)
        timeSinceSwitch = 10;
        uInd_history = 1;
    else
        timeSinceSwitch = params.simu.time(timeInd) - params.simu.time(lastSwitchInd);
    end


% Get position and velocity
thetaDot = states(1);
theta = states(2);

% Get torque options
ptoTorqueOptions = params.hyd.Force2Torque(theta)*params.hyd.ptoForceOptions(:);
if ctrl.limitChoices
    % ptoTorqueOptions([5,6,8,9]) = NaN;
end


% Sliding Mode Contrl
thetaError = theta-ctrl.optTraj.theta(timeInd);
thetaDotError = thetaDot-ctrl.optTraj.thetaDot(timeInd);
s = thetaDotError+ctrl.lambda*thetaError;

% Precompute Transition Matrices
[M,H] = getTransition(params,2*ctrl.horizonInd);

% pad vectors so we can simulate past the official time vector
thetaOpt = [ctrl.optTraj.theta;zeros(2*ctrl.horizonInd,1)];
thetaDotOpt = [ctrl.optTraj.thetaDot;zeros(2*ctrl.horizonInd,1)];
Texc = [wave.torque.Texc;zeros(2*ctrl.horizonInd,1)];

% if s is large, recalulate control input
if abs(s) > ctrl.phi || timeInd == 1
    % Horizon indices
    hInds = (timeInd:(timeInd+2*ctrl.horizonInd-1));

    % Desired States (we dont care about what happens with the rad states
    xStar = [thetaDotOpt(hInds),thetaOpt(hInds),zeros(length(hInds),2)];

    % Initialize cost vector
        % Decaying function from previous switching time
    J = (exp(-10*timeSinceSwitch) + ctrl.phi) * ones(length(ptoTorqueOptions),length(ptoTorqueOptions));
    
        % Constant loss if switching at all
    J = J .* ~eye(size(J));

        % Zero switching loss if we start from the choice we are currently at
    J(uInd_history(end),:) = 0;

        

    % Simulate future times for each control input
    for uInd = 1:length(ptoTorqueOptions)
        % Control inputs:
        uOption = ptoTorqueOptions(uInd) * ones(ctrl.horizonInd,1);
        for uInd2 = 1:length(ptoTorqueOptions)
        % Control inputs:
        uOption(ctrl.horizonInd+1:2*ctrl.horizonInd) = ptoTorqueOptions(uInd2) * ones(ctrl.horizonInd,1);
        
        % simulate this control action over the time horizon
        xHat_stacked = M*states + H*(uOption+Texc(hInds));
        xHat = reshape(xHat_stacked,size(xStar,2),size(xStar,1))';

        % xHat2 = lsim(params.phys.sys,uOption+Texc(hInds),params.simu.time(hInds));

        S = (xHat - xStar)*[1; ctrl.lambda; 0; 0];

        % J(uInd) = J(uInd) + sum(S.^2);
        J(uInd,uInd2) = J(uInd,uInd2) + sum(S.^2);
        % J(uInd,uInd2) = J(uInd,uInd2) + abs(S(end));
        end
    end
    [~,tmp] = min(J,[],"all");
    [uOptInd,~] = ind2sub(size(J),tmp);
else
    uOptInd = uInd_history(end);
end

% Output 
out.controlValue = ptoTorqueOptions(uOptInd);
out.controlIndex = uOptInd;
end

function [M,H] = getTransition(params,N)

Phi = eye(size(params.phys.sys.A)) + params.phys.sys.A * params.simu.dt;
Gamma = params.phys.sys.B * params.simu.dt;

nx = size(Phi, 1);

% Pre-allocate the large matrices
M = zeros(nx * N, nx);
H = zeros(nx * N, N);

% Precompute the first column of H and the rows of M
% We iterate forward: Phi^1*Gamma, Phi^2*Gamma, etc.
current_Phi = Phi; 
current_Gamma = Gamma;

for i = 1:N
    row_idx = (i-1)*nx + 1 : i*nx;
    
    % Fill M (Initial state effect)
    M(row_idx, :) = current_Phi;
    
    % Fill the first column-block of H_big
    H(row_idx, 1) = current_Gamma;
    
    % Update for next step
    current_Phi = Phi * current_Phi;
    current_Gamma = Phi * current_Gamma;
end

% Fill the rest of H
% Each column is just a shifted version of the column to its left
for j = 2:N
    % Copy the previous column shifted down by nx rows
    source_rows = 1 : (N-j+1)*nx;
    dest_rows = (j-1)*nx + 1 : N*nx;
    
    source_cols = (j-2) + 1 : (j-1);
    dest_cols = (j-1) + 1 : j;

    source_cols = 1;
    dest_cols = j;
    
    H(dest_rows, dest_cols) = H(source_rows, source_cols);
end
end