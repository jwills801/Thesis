function ctrl = getControl(params,wave)



% controller = 'PI';
% controller = 'Sliding Mode'
% controller = 'Coulomb Damping';
%controller = 'MPC_QP';
controller = 'MPC_DP';
ctrl.controller = controller;

% Define parameters for all controllers
ctrl.timeHorizon = .2; % Length of a control step
ctrl.horizonInd = round(ctrl.timeHorizon/params.simu.dt);
ctrl.numHorizons = 100;

% Get optimal trajectory and energy
optTraj = getOptimal(params,wave);
ctrl.optTraj = optTraj;

% Define parameters unique to each controller
switch controller
    case 'PI'
        ctrl.limitChoices = 1;
    case 'Sliding Mode'
        ctrl.lambda = 1; % Defines the sliding surface
        ctrl.phi = 3e-2; % band around sliding surface
    case 'Coulomb Damping'
    case'MPC_QP'
        m = ctrl.numHorizons;
        gamma = 3e-6;

        % Precompute Transition Matrices
        [M,H] = getTransition(params,m*ctrl.horizonInd);

        % Other Matrices
        [L,C,Q_sw] = getUtilityMatrices(m,ctrl.horizonInd);

        Q = transpose(H*L)*C + gamma*Q_sw;

        % output
        ctrl.MPC.m = m;
        ctrl.MPC.gamma = gamma;
        ctrl.MPC.M = M;
        ctrl.MPC.H = H;
        ctrl.MPC.L = L;
        ctrl.MPC.C = C;
        ctrl.MPC.Q = Q;

    case'MPC_DP'
        m = ctrl.numHorizons;
        n = ctrl.horizonInd;
        N = m*n;
        dt= params.simu.dt;
        m_Astar = 5;

        % Precompute Matrices for mechanical energy calculation
        [M_local,H_local] = getTransition(params,n);
        [~,C_local,~] = getUtilityMatrices(1,n);
        ctrl.w = M_local'*C_local*dt;
        ctrl.Q_local = ones(1,n)*H_local'*C_local*dt;
        ctrl.b_exc = waveEnergyContribution(C_local'*H_local*dt,wave.torque.Texc);

        % Precompute terminal cost matrices
        [M,H] = getTransition(params,N);
        [L,C,~] = getUtilityMatrices(m,n);
        Hu = H*L;
        Qinv = inv(Hu'*C + C'*Hu);
        ctrl.P=1/2*dt*M'*C*Qinv*C'*M;
        [ctrl.rT,ctrl.b] = waveTerminalCost(H,C,Qinv,M,dt,wave.torque.Texc);

        % Compute uncontroled trajectory
        ctrl.Xfree = lsim(params.phys.sys,wave.torque.Texc,wave.torque.time);
end

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

function [L,C,Q_sw] = getUtilityMatrices(m,n)
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

% Make Q_sw
main_diag = 2 * ones(m, 1);
main_diag(end) = 1;

off_diag = -1 * ones(m-1, 1);

% diag(v, k) places vector v on the k-th diagonal
Q_sw = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
end

function b_exc = waveEnergyContribution(gain,Texc)
N = size(Texc,1); % Number of fine time steps
n = size(gain,2); % Number of fine time steps in one coarse one
m = floor(N/n); % Number of coarse time steps

b_exc = NaN(m,1);
for k = 1:m
    % Extract the n-step window for this specific block k
    idx_start = (k-1)*n + 1;
    idx_end   = k*n;

    % T_k contains the n wave torque samples for the current block
    T_k = Texc(idx_start:idx_end);

    % Calculate scalar wave energy potential for this block
    b_exc(k) = gain * T_k;
end
end

function [rT, b] = waveTerminalCost(H,C,Qinv,M,dt,Texc)
% Recover number of control time steps
N = size(H,2); % Number of fine time steps to be considered in the terminal cost
m = size(Qinv,1); % Number of coarse time steps in the terminal cost
n = floor(N/m); % Number of fine time steps in one coarse one

% Total number of fine and coarse time steps
N_total = size(Texc,1); % Number of fine time steps
m_total = floor(N_total/n)-m;

% Group all terms that don't depend on the specific wave window T
K_r = H' * C * Qinv * C' * M;     % [N x 4] matrix
K_b = H' * C * Qinv * C' * H; % [N x N] matrix

% Initilize vectors
rT = NaN(m_total,4);
b = NaN(m_total,1);

% Calculate rT and b for each control window
for k = 1:m_total
    % Extract the N-step window for this specific block k
    idx_start = (k-1)*n + 1;
    idx_end   = idx_start + N - 1;
    
    % T contains the N wave torque samples for the current block
    T = Texc(idx_start:idx_end);
    
    % Calculate rT and b
    rT(k,:) = T' * K_r * dt;   % Result is [1 x 4]
    b(k)    = T' * K_b * T * dt/2; % Result is scalar
end

end