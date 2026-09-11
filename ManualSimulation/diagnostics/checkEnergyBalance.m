% checkEnergyBalance.m
% Runs a short simulation for each drivetrain family and checks that
% average mechanical power in equals average electrical power out plus
% average total losses (aveMechPow == aveElecPow + aveLoss). Before the
% evaluate.m fix this fails for DHD/PassivePump by ~2.25% of aveElecPow;
% after the fix it should pass for all three families.
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'wave'));
addpath(fullfile(root,'control')); addpath(fullfile(root,'dynamics'));
addpath(fullfile(root,'evaluation'));

cases = {
    'PassivePump','CoulombDamping',2
    'EHA',        'MPC_QP',        0
    'DHD',        'MPC_Astar',     2
};

fprintf('--- checkEnergyBalance ---\n');
for i = 1:size(cases,1)
    runParams = struct();
    runParams.drive = cases{i,1};
    runParams.controller = cases{i,2};
    runParams.pressure_rails = cases{i,3};
    switch runParams.drive
        case 'PassivePump'
            runParams.rodArea = (0.0254*6)^2*pi;
            runParams.capArea = 1.5*runParams.rodArea;
            runParams.highPressure = 20.6e6;
        case 'EHA'
            runParams.considerLosses = 0;
            runParams.rodArea = (0.0254*8)^2*pi;
            runParams.capArea = runParams.rodArea;
        case 'DHD'
            runParams.considerLosses = 1;
            runParams.rodArea = (.0254*6)^2*pi;
            runParams.capArea = 1.5*runParams.rodArea;
            runParams.highPressure = 33.3e6;
    end

    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.finalTime = 60; % short run for diagnostics speed
    params.simu.rampTime = 10;
    params.simu.time = (0:params.simu.dt:params.simu.finalTime)';

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    err = abs(ev.aveMechPow - (ev.aveElecPow + ev.aveLoss));
    fprintf('%s / %s: ', runParams.drive, runParams.controller);
    report('aveMechPow == aveElecPow + aveLoss', err, max(1,abs(ev.aveMechPow))*1e-6);
end

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
