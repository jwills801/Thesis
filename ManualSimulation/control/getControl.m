function ctrl = getControl(params,wave)

%% Define parameters for all controllers
ctrl.timeHorizon = .2; % Length of a control step
ctrl.horizonInd = round(ctrl.timeHorizon/params.simu.dt);
ctrl.numHorizons = params.simu.peakPeriod*2.5/ctrl.timeHorizon;

%% Get optimal trajectory and energy
optTraj = getOptimal(params,wave);
ctrl.optTraj = optTraj;

%% Define parameters unique to each controller
switch params.runParams.controller
    case 'Coulomb Damping'
    case 'PI'
    case'MPC_QP'
        % Precompute Transition Matrices
        m = ctrl.numHorizons;
        [M,H] = getTransition(params,m*ctrl.horizonInd);

        % Other Matrices
        [L,C] = getUtilityMatrices(m,ctrl.horizonInd);

        % save output to a structure
        ctrl.MPC.m = m;
        ctrl.MPC.M = M;
        ctrl.MPC.H = H;
        ctrl.MPC.L = L;
        ctrl.MPC.C = C;

        switch params.runParams.drive
            case 'EHA'
                ctrl.MPC.EHA = MPC_EHA(params,ctrl.MPC);
        end

    case 'MPC_Astar'
        % Number of A star time steps
        ctrl.m_Astar = 5;

        % unwrap useful parameters
        m = ctrl.numHorizons; % This is the terminal cost horizon
        n = ctrl.horizonInd;
        dt = ctrl.timeHorizon;
        
        % Precompute Matrices for mechanical energy calculation
        [M_local,H_local] = getTransition(params,n);
        [~,C_local] = getUtilityMatrices(1,n);
        ctrl.w = M_local'*C_local*dt;
        ctrl.Q_local = ones(1,n)*H_local'*C_local*dt;
        ctrl.b_exc = waveEnergyContribution(C_local'*H_local*dt,wave.torque.Texc);
        T_block_full = waveEnergyContribution(H_local,wave.torque.Texc);
        H_block_full = H_local*ones(n,1);
        
        % Just take the last 4 rows
        ctrl.T_block = T_block_full(end-3:end,:);
        ctrl.H_block = H_block_full(end-3:end);
        ctrl.M_block = M_local(end-3:end,:);

        % Precompute terminal cost matrices
        ctrl.termCost = getTerminalCostMatrices(params,wave,dt,m,n,ctrl.m_Astar);

        % Test the matrices
        % testing(ctrl,params,wave);
end

end




function [M,H] = getTransition(params,N)

% Forward Euler
Phi = eye(size(params.phys.sys.A)) + params.phys.sys.A * params.simu.dt;
Gamma = params.phys.sys.B * params.simu.dt;

% Zero Order Hold
A = params.phys.sys.A;
B = params.phys.sys.B;
dt = params.simu.dt;
Phi = expm(A * dt);
Gamma = integral(@(s) expm(A*s), 0, dt, 'ArrayValued', true) * B;

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
N_total = size(Texc,1); % Total number of fine time steps
m_total = floor(N_total/n)-m; % total number of coarse time steps

% Group all terms that don't depend on the specific wave window T
K_r = H' * C * Qinv * C' * M; % [N x 4] matrix
K_b = H' * C * Qinv * C' * H; % [N x N] matrix

% Initilize vectors
rT = NaN(m_total,4);
b = NaN(m_total,1);

% Calculate rT and b for each control window
for k = 1:m_total
    % Extract the N-step window for this specific block k
    idx_start = (k-1)*n + 1; % if k = 1, start at zero
    idx_end   = idx_start + N - 1;
    
    % T contains the N wave torque samples for the current block
    T = Texc(idx_start:idx_end);
    
    % Calculate rT and b
    rT(k,:) = -T' * K_r * dt;   % Result is [1 x 4]
    b(k)    = -T' * K_b * T * dt/2; % Result is scalar
end

end

function out = getTerminalCostMatrices(params,wave,dt,m,n,m_Astar)
% Comput terminal cost matrices for inside the Astar method
out(m_Astar) = struct();
bar = waitbar(0,'Precomputing Terminal Cost Matrices');
for d = 1:m_Astar
    waitbar((d-1)/m_Astar,bar)
    N = (m-d)*n;

    % Get matrices
    [M,H] = getTransition(params,N);
    [L,C] = getUtilityMatrices(m-d,n);
    Hu = H*L;
    Qinv = inv(Hu'*C + C'*Hu);
    P = -1/2*dt*M'*C*Qinv*C'*M;
    [rT,b] = waveTerminalCost(H,C,Qinv,M,dt,wave.torque.Texc);

    % Save matrices
    out(d).P=P;
    out(d).rT = rT;
    out(d).b = b;
end
close(bar)
end

function out = MPC_EHA(params,in)

% Choose whether to include losses in the controller
if params.runParams.considerLosses
    a = params.hyd.EHA.LossCoeffs(1);
    b = params.hyd.EHA.LossCoeffs(2);
    c = params.hyd.EHA.LossCoeffs(3);
    d = params.hyd.EHA.LossCoeffs(4);
else
    a = 0; b = 0; c = 0; d = 0;
end

% Construct matrices from the loss coeffs
out = struct();
out.Q = in.C * diag(a*ones(in.m,1)) * in.C';
out.G = in.C * diag(b*ones(in.m,1));
out.R = diag(c*ones(in.m,1));
out.D = d;

end


function testing(ctrl,params,wave)

%% simulate trajectory
t_start = 50; [~,t_startInd] = min(abs(wave.torque.time-t_start));
t_end = t_start + 0.2; [~,t_endInd] = min(abs(wave.torque.time-t_end));
time = wave.torque.time(t_startInd:t_endInd);
Texc = wave.torque.Texc(t_startInd:t_endInd);
x0 = zeros(4,1);
sys = params.phys.sys;

% Compute optimal control
k=floor(t_start/.2)+1;
w = ctrl.w;
b_exc = ctrl.b_exc(k);
Q_local = ctrl.Q_local;
u = -(w'*x0+b_exc)/2/Q_local;
% u = 1e6;%*sin(time);
u = 0;

% Advance state
M_block = ctrl.M_block;
T_block = ctrl.T_block(:,k);
H_block = ctrl.H_block;
xf = M_block*x0+T_block+H_block*u;

% Energy in this block
E_mech = (w'*x0+b_exc)*u + u^2*Q_local;

% future cost
d = 1;
P = ctrl.termCost(d).P;
rT = ctrl.termCost(d).rT(k+1,:);
b = ctrl.termCost(d).b(k+1);
E_term1 = xf'*P*xf + rT*xf + b
xf'*P*xf

% Now lets look at the next time step
x0 = xf; d=2;
% choose optimal control again
b_exc = ctrl.b_exc(k+d-1);
u = -(w'*x0+b_exc)/2/Q_local;
u = 0*16.096e6;
% Advance state again
T_block = ctrl.T_block(:,k+d-1);
xf = M_block*x0+T_block+H_block*u;
% Compute absorbed energy
E_mech = (w'*x0+b_exc)*u + u^2*Q_local
% compute terminal energy
rT = ctrl.termCost(d).rT(k+d,:);
b = ctrl.termCost(d).b(k+d);
E_term = xf'*P*xf + rT*xf + b
xf'*P*xf
% total cost
E_total = E_mech+E_term
E_total-E_term1

disp('---------------')
%% Compare these states to lsim
% Lsim results
states = lsim(sys,Texc+u,time,x0, 'zoh');
mechPow = u.*states(2:end,1);
mechEnergy = sum(mechPow)*(time(2)-time(1));

% Compare results
[xf(1:2),states(end,1:2)'];
[xf(3:4),states(end,3:4)'];

disp('---------------')
% figure, plot(time,states(:,1))
a=0;
end