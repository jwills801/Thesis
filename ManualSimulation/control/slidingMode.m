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

% Sliding Mode Contrl
thetaError = theta-ctrl.optTraj.theta(timeInd);
thetaDotError = thetaDot-ctrl.optTraj.thetaDot(timeInd);
s = thetaDotError+ctrl.lambda*thetaError;

% Precompute Transition Matrices
[M,H] = getTransition(params,ctrl.horizonInd);

% pad vectors so we can simulate past the official time vector
thetaOpt = [ctrl.optTraj.theta;zeros(ctrl.horizonInd,1)];
thetaDotOpt = [ctrl.optTraj.thetaDot;zeros(ctrl.horizonInd,1)];
Texc = [wave.torque.Texc;zeros(ctrl.horizonInd,1)];

% if s is large, recalulate control input
if abs(s) > ctrl.phi || timeInd == 1
    % Horizon indices
    hInds = (timeInd:(timeInd+ctrl.horizonInd-1));

    % Desired States (we dont care about what happens with the rad states
    xStar = [thetaDotOpt(hInds),thetaOpt(hInds),zeros(length(hInds),2)];

    % Initialize cost vector
    J = exp(-10*timeSinceSwitch) * ones(size(ptoTorqueOptions));
    J(uInd_history(end)) = 0;

    % Simulate future times for each control input
    for uInd = 1:length(ptoTorqueOptions)
        % Control inputs:
        uOption = ptoTorqueOptions(uInd);
        
        % simulate this control action over the time horizon
        xHat_stacked = M*states + H*(uOption+Texc(hInds));
        xHat = reshape(xHat_stacked,size(xStar,2),size(xStar,1))';

        S = (xHat - xStar)*[1; ctrl.lambda; 0; 0];

        % J(uInd) = J(uInd) + sum(S.^2);
        J(uInd) = J(uInd) + abs(S(end));
    end
    [~,uOptInd] = min(J);
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
    
    H(dest_rows, dest_cols) = H(source_rows, source_cols);
end
end