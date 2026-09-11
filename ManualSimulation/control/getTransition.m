% getTransition.m
% Builds the block-Toeplitz zero-order-hold prediction matrices M,H from
% params.phys.sys: x[1:N] = M*x0 + H*u[0:N-1].
% Note: control/slidingMode.m has its OWN local function of the same name
% (a separate forward-Euler implementation) that shadows but does NOT
% call this file -- not a real dependency edge, just a naming collision.
% Calls: none
% Called by: control/getControl.m, diagnostics/checkMPC_EHA.m,
%   checkTerminalCost.m, checkTransitionMatrices.m
function [M,H] = getTransition(params,N)
% Builds the block transition matrices for the zero-order-hold discretized
% WEC state space, unrolled over N fine time steps:
%   x[1:N] = M*x[0] + H*(u[0:N-1])
% M stacks Phi^1..Phi^N (initial-state propagation).
% H is the block-Toeplitz input-to-state map built from Phi^k*Gamma.
A = params.phys.sys.A;
B = params.phys.sys.B;
dt = params.simu.dt;

Phi = expm(A * dt);
Gamma = integral(@(s) expm(A*s), 0, dt, 'ArrayValued', true) * B;

nx = size(Phi, 1);

M = zeros(nx * N, nx);
H = zeros(nx * N, N);

current_Phi = Phi;
current_Gamma = Gamma;

for i = 1:N
    row_idx = (i-1)*nx + 1 : i*nx;
    M(row_idx, :) = current_Phi;
    H(row_idx, 1) = current_Gamma;

    current_Phi = Phi * current_Phi;
    current_Gamma = Phi * current_Gamma;
end

% Each later column of H is the first column shifted down by nx rows
for j = 2:N
    source_rows = 1 : (N-j+1)*nx;
    dest_rows = (j-1)*nx + 1 : N*nx;
    H(dest_rows, j) = H(source_rows, 1);
end
end
