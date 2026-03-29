function ctrl = getControl(params,wave)

% controller = 'PI';
    ctrl.limitChoices = 1;
% controller = 'Sliding Mode';
    ctrl.timeHorizon = .2; % Amount of time to simulate
    ctrl.horizonInd = round(ctrl.timeHorizon/params.simu.dt);
    ctrl.lambda = 1; % Defines the sliding surface
    ctrl.phi = 3e-2; % band around sliding surface
% controller = 'Coulomb Damping';
controller = 'MPC';

optTraj = getOptimal(params,wave);

ctrl.optTraj = optTraj;
ctrl.controller = controller;

switch controller
    case'MPC'
    % numberof control time steps
m = 100;
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