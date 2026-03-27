function out = MPC(params,ctrl,wave,states,uInd_history)
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

m = 10;

% Precompute Transition Matrices
[M,H] = getTransition(params,m*ctrl.horizonInd);

% Other Matrices
[L,C] = getUtilityMatrices(m,ctrl.horizonInd);

Q = transpose(H*L)*C;

% pad vectors so we can simulate past the official time vector
Texc = [wave.torque.Texc;zeros(m*ctrl.horizonInd,1)];

% place holder u values
u = ptoTorqueOptions(ones(m,1));

% switch on timehorizon time intervals
if 1% mod(timeInd,ctrl.horizonInd) == 1
    % Horizon indices
    hInds = (timeInd:(timeInd+m*ctrl.horizonInd-1));

    X = M*states + H*Texc(hInds) + H*L*u;

    E = transpose(X)*C*u*.01;

    u = inv(-Q+transpose(-Q)) * transpose(C)*(M*states+H*Texc(hInds));
else
    u = uInd_history(end);
end
% Output 
out.controlValue = u(1);
out.controlIndex = 1;
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

function [L,C] = getUtilityMatrices(m,n)
    onesCol = ones(n,1);
    zerosCol = zeros(n,1);
    
    ei = [1;0;0;0];
    c_block = repmat(ei,n,1);

    L = NaN(n*m,m);
    C = zeros(4*n*m,m);
    for col = 1:m
        L(:,col) = [repmat(zerosCol,col-1,1);
                    onesCol;
                    repmat(zerosCol,m-col,1)];
        C(:,col) = [repmat(0*c_block,col-1,1);
                    c_block;
                    repmat(0*c_block,m-col,1)];
    end
    
end