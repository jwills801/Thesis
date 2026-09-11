% checkTerminalCost.m
% Verifies ctrl.termCost(d).{P,rT,b} (the A* search's terminal-cost
% heuristic, from getControl.m's getTerminalCostMatrices/waveTerminalCost)
% against brute force: propagate the given state forward N fine steps
% under the closed-form optimal continuous control u*=-Q^{-1}f and
% directly compute the resulting absorbed energy X'Cu*dt.
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, getTransition.m, getUtilityMatrices.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'wave'));
addpath(fullfile(root,'control'));

runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',2, ...
    'considerLosses',1,'rodArea',(0.0254*6)^2*pi,'capArea',1.5*(0.0254*6)^2*pi, ...
    'highPressure',35e6);
params = getParameters(runParams);
% NOTE: ctrl.numHorizons = 2.5*peakPeriod/timeHorizon (timeHorizon=0.2 for
% MPC_Astar) is never rounded in getControl.m -- a non-integer value
% breaks NaN(n*m,m) downstream in getUtilityMatrices. Found via this
% diagnostic (see diagnostics/ReadMe.md). Picking peakPeriod=0.64 keeps
% numHorizons=8 exactly integer (with room for m_Astar=5 terminal-cost
% horizons, m-d>=1) while still shrinking the matrices for a fast check.
params.simu.peakPeriod = 0.64;
wave = generateExcitingTorque(params);
ctrl = getControl(params,wave);

d = 1; % first terminal-cost horizon
m = ctrl.numHorizons; n = ctrl.horizonInd;
N = (m-d)*n;
[M,H] = getTransition(params,N);
[L,C] = getUtilityMatrices(m-d,n);
Hu = H*L;
Qinv = inv(Hu'*C + C'*Hu);

x0 = [0.05;-0.02;0.01;-0.005];
k = 3; % arbitrary control block index within range
T = wave.torque.Texc((k-1)*n+1 : (k-1)*n+N);
dt = params.simu.dt;

f = C'*M*x0 + C'*H*T;
uStar = -Qinv*f;
X = M*x0 + H*T + Hu*uStar;
E_bruteforce = (X'*C*uStar)*dt;

P = ctrl.termCost(d).P;
rT = ctrl.termCost(d).rT(k,:);
b = ctrl.termCost(d).b(k);
E_closedform = x0'*P*x0 + rT*x0 + b;

err = abs(E_bruteforce - E_closedform);
fprintf('--- checkTerminalCost ---\n');
report('Closed-form terminal cost matches brute-force optimal-control energy', ...
    err, max(1,abs(E_bruteforce))*1e-6);

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
