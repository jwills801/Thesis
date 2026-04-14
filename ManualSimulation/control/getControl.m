function ctrl = getControl(params,wave)



% controller = 'PI';
% controller = 'Sliding Mode'
% controller = 'Coulomb Damping';
% controller = 'MPC_QP';
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
        gamma = 0*3e-6;

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
        dt= params.simu.dt;
        ctrl.m_Astar = 5;

        % Precompute Matrices for mechanical energy calculation
        [M_local,H_local] = getTransition(params,n);
        [~,C_local,~] = getUtilityMatrices(1,n);
        ctrl.w = M_local'*C_local*dt;
        ctrl.Q_local = ones(1,n)*H_local'*C_local*dt;
        ctrl.b_exc = waveEnergyContribution(C_local'*H_local*dt,wave.torque.Texc);
        ctrl.Q_local2 = ctrl.Q_local + ctrl.Q_local';

        T_block_full = waveEnergyContribution(H_local,wave.torque.Texc);
        H_block_full = H_local*ones(n,1);
        
        % Just take the last 4 rows
        ctrl.T_block = T_block_full(end-3:end,:);
        ctrl.H_block = H_block_full(end-3:end);
        ctrl.M_block = M_local(end-3:end,:);

        % Precompute terminal cost matrices
        ctrl.termCost = getTerminalCostMatrices(params,wave,m,n,ctrl.m_Astar);

        % Test matrices
       %  testing(ctrl,params,wave);
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

function out = waveEnergyContribution(gain,Texc)
N = size(Texc,1); % Number of fine time steps
n = size(gain,2); % Number of fine time steps in one coarse one
m = floor(N/n); % Number of coarse time steps

% This function works for both b_exc and T_block
out = NaN(size(gain,1),m);
for k = 1:m
    % Extract the n-step window for this specific block k
    idx_start = (k-1)*n + 1;
    idx_end   = k*n;

    % T_k contains the n wave torque samples for the current block
    T_k = Texc(idx_start:idx_end);

    % Calculate scalar wave energy potential for this block
    out(:,k) = gain * T_k;
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
    rT(k,:) = -T' * K_r * dt;   % Result is [1 x 4]
    b(k)    = -T' * K_b * T * dt/2; % Result is scalar
end

end

function out = getTerminalCostMatrices(params,wave,m,n,m_Astar)
% Comput terminal cost matrices for inside the Astar method
dt= params.simu.dt;
out(m_Astar) = struct();
for d = 1:m_Astar
    N = (m-d)*n;

    % Get matrices
    [M,H] = getTransition(params,N);
    [L,C,~] = getUtilityMatrices(m-d,n);
    Hu = H*L;
    Qinv = inv(Hu'*C + C'*Hu);
    P = -1/2*dt*M'*C*Qinv*C'*M;
    [rT,b] = waveTerminalCost(H,C,Qinv,M,dt,wave.torque.Texc);

    % Save matrices
    out(d).P=P;
    out(d).rT = rT;
    out(d).b = b;
end
end


function testing(ctrl,params,wave)

        %% simulate trajectory
            t_start = 50; [~,t_startInd] = min(abs(wave.torque.time-t_start));
            t_end = t_start + 0.2; [~,t_endInd] = min(abs(wave.torque.time-t_end));
            time = wave.torque.time(t_startInd:t_endInd);
            Texc = wave.torque.Texc(t_startInd:t_endInd);
            x0 = zeros(4,1);
            sys = params.phys.sys;
            u = 1e6;%*sin(time);
        
            states = lsim(sys,Texc+u,time,x0);

            mechPow = u.*states(:,1);
            mechEnergy = trapz(time,mechPow)

            % Test matrices
                    % Advance state
                    k=floor(t_start/.2)+1;
        M_block = ctrl.M_block;
        T_block = ctrl.T_block(:,k);
        H_block = ctrl.H_block;
        xf = M_block*x0+T_block+H_block*u;
        [xf(1:2),states(end,1:2)']

        % Energy in this block
        w = ctrl.w;
        b_exc = ctrl.b_exc(k);
        Q_local = ctrl.Q_local;
        E_mech = (w'*x0+b_exc)*u + u^2*Q_local

        % future cost
        P = ctrl.termCost(1).P;
        rT = ctrl.termCost(1).rT(k,:);
        b = ctrl.termCost(1).b(k);
        E_term = xf'*P*xf + rT*xf + b


% figure, plot(time,states(:,1))
a=0;
end