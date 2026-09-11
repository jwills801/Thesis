% checkMPC_EHA.m
% Cross-checks getControl's ctrl.MPC.{A,Bx,Bt} (used by MPC_QP.m to solve
% u* = -(A+A')\B', with B = x0'*Bx + T'*Bt) against a finite-difference
% gradient of the electrical-energy cost, built independently here from
% the main(1).tex EHA-loss/"Continuous Quadratic Programming" derivation
% (E_elec = X'Cu + X'QX + X'GU + U'RU + D), not by calling MPC_EHA itself.
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, getTransition.m, getUtilityMatrices.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'wave'));
addpath(fullfile(root,'control'));

for considerLosses = [0 1]
    runParams = struct('drive','EHA','controller','MPC_QP', ...
        'considerLosses',considerLosses, ...
        'rodArea',(0.0254*8)^2*pi,'capArea',(0.0254*8)^2*pi);
    params = getParameters(runParams);
    params.simu.peakPeriod = 0.2; % shrink horizon size for a fast, small-scale check
    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);

    m = ctrl.numHorizons; n = ctrl.horizonInd; N = m*n;
    [M,H] = getTransition(params,N);
    [L,C] = getUtilityMatrices(m,n);
    Hu = H*L;

    a=0; b=0; c=0; d=0; %#ok<NASGU>
    if considerLosses
        a = params.hyd.EHA.LossCoeffs(1); b = params.hyd.EHA.LossCoeffs(2);
        c = params.hyd.EHA.LossCoeffs(3); d = params.hyd.EHA.LossCoeffs(4); %#ok<NASGU>
    end
    Q = diag(repmat([a 0 0 0],1,N));
    G = C*diag(b*ones(m,1));
    R = L'*diag(c*ones(N,1))*L;

    x0 = [0.05;-0.02;0.01;-0.005];
    T = wave.torque.Texc(1:N);

    Efun = @(u) energyElec(M,H,Hu,C,Q,G,R,x0,T,u);

    u0 = 1e5*(1:m)'; % arbitrary nonzero point, away from any optimum
    grad_numeric = numericGradient(Efun,u0);

    B = x0'*ctrl.MPC.Bx + T'*ctrl.MPC.Bt;
    grad_analytic = (ctrl.MPC.A+ctrl.MPC.A')*u0 + B';

    err = max(abs(grad_numeric(:)-grad_analytic(:)));
    fprintf('--- checkMPC_EHA (considerLosses=%d) ---\n',considerLosses);
    report('ctrl.MPC gradient matches finite-difference of E_elec(u)', ...
        err, max(1,norm(grad_analytic))*1e-4);
end

function E = energyElec(M,H,Hu,C,Q,G,R,x0,T,u)
X = M*x0 + H*T + Hu*u;
E = X'*C*u + X'*Q*X + X'*G*u + u'*R*u;
end

function g = numericGradient(f,u0)
h = 1;
g = zeros(size(u0));
for i = 1:length(u0)
    up = u0; up(i) = up(i)+h;
    um = u0; um(i) = um(i)-h;
    g(i) = (f(up)-f(um))/(2*h);
end
end

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
