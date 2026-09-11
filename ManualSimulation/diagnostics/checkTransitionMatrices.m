% checkTransitionMatrices.m
% Verifies control/getTransition.m's M,H against brute-force repeated
% zero-order-hold propagation, using distinct nonzero inputs at each step
% so a column-ordering bug would show up.
% Calls: parameters/getParameters.m, control/getTransition.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'control'));

runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',2, ...
    'considerLosses',1,'rodArea',(0.0254*6)^2*pi,'capArea',1.5*(0.0254*6)^2*pi, ...
    'highPressure',35e6);
params = getParameters(runParams);

N = 6; % small horizon for brute-force comparison
[M,H] = getTransition(params,N);

A = params.phys.sys.A; B = params.phys.sys.B; dt = params.simu.dt;
Phi = expm(A*dt);
Gamma = integral(@(s) expm(A*s),0,dt,'ArrayValued',true)*B;

nx = size(A,1);
x0 = [0.1;-0.05;0.02;-0.01]; % arbitrary nonzero IC
u = (1:N)'*1e5;              % arbitrary distinct inputs

xBrute = zeros(nx,N);
x = x0;
for k = 1:N
    x = Phi*x + Gamma*u(k);
    xBrute(:,k) = x;
end

xFromMH = reshape(M*x0 + H*u, nx, N);

err = max(abs(xFromMH(:) - xBrute(:)));
fprintf('--- checkTransitionMatrices ---\n');
report('M,H reproduce brute-force ZOH propagation',err,1e-6);

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g)\n',name,err);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
