% validatePhase2Subset.m
% End-to-end validation of the phase-2 machinery (sea-state loop +
% per-sea-state pressure grid search) on a small subset: 2 placeholder
% Humboldt sea-state bins x a coarse pressure grid x one case per
% drivetrain family, with a short finalTime to keep runtime down. Not the
% production sweep -- just proves the harness produces sane, non-NaN
% results before a full run.
% Calls: wave/humboldtSeaStates.m, parameters/optimizePressure.m,
%   getParameters.m, generateExcitingTorque.m, control/getControl.m,
%   dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters')); addpath(fullfile(root,'wave'));
addpath(fullfile(root,'control')); addpath(fullfile(root,'dynamics'));
addpath(fullfile(root,'evaluation'));

seaStates = humboldtSeaStates();
seaStates = seaStates(1:2,:); % validation subset: first 2 bins only

simOverrides = struct('finalTime',60,'rampTime',10);
simOverrides.time = (0:0.01:simOverrides.finalTime)';

fprintf('=== validatePhase2Subset ===\n');
for s = 1:height(seaStates)
    seaState = struct('Hs',seaStates.Hs(s),'Tp',seaStates.Tp(s));
    fprintf('\n--- Sea state %d: Hs=%.2gm, Tp=%.2gs ---\n',s,seaState.Hs,seaState.Tp);

    % PassivePump / CoulombDamping (no switchMap regen needed)
    runParams = struct('drive','PassivePump','controller','CoulombDamping', ...
        'pressure_rails',2,'rodArea',(0.0254*6)^2*pi,'capArea',1.5*(0.0254*6)^2*pi);
    res = optimizePressure(runParams,seaState,[20e6 35e6],simOverrides);
    reportSweep('PassivePump/CoulombDamping',res);

    % EHA / MPC_QP, both mechanical- and electrical-optimized (no pressure rails)
    for considerLosses = [0 1]
        thisRunParams = struct('drive','EHA','controller','MPC_QP', ...
            'considerLosses',considerLosses,'rodArea',(0.0254*8)^2*pi,'capArea',(0.0254*8)^2*pi);
        params = getParameters(thisRunParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = seaState.Hs;
        params.simu.peakPeriod = seaState.Tp;
        params.simu.finalTime = simOverrides.finalTime;
        params.simu.rampTime = simOverrides.rampTime;
        params.simu.time = simOverrides.time;

        wave = generateExcitingTorque(params);
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);
        ev = evaluate(params,dyn,ctrl);
        fprintf('EHA/MPC_QP (considerLosses=%d): elecRGP=%.3g mechRGP=%.3g %s\n', ...
            considerLosses,ev.elecRGP,ev.mechRGP,validity(ev.elecRGP));
    end

    % DHD / MPC_Astar, 2 rails (WITH switchMap regen per grid point -- slow)
    runParams = struct('drive','DHD','controller','MPC_Astar', ...
        'pressure_rails',2,'considerLosses',1, ...
        'rodArea',(0.0254*6)^2*pi,'capArea',1.5*(0.0254*6)^2*pi);
    res = optimizePressure(runParams,seaState,[30e6 35e6],simOverrides);
    reportSweep('DHD/MPC_Astar (2 rails)',res);
end

function s = validity(x)
if isnan(x) || ~isfinite(x)
    s = '<-- FAIL (NaN/Inf)';
else
    s = '';
end
end

function reportSweep(name,res)
for i = 1:numel(res.pressureGrid)
    fprintf('%s @ %.3gMPa: elecRGP=%.3g mechRGP=%.3g %s\n', ...
        name,res.pressureGrid(i)/1e6,res.elecRGP(i),res.mechRGP(i),validity(res.elecRGP(i)));
end
fprintf('%s BEST: %.3gMPa (elecRGP=%.3g)\n',name,res.bestPressure/1e6,res.bestElecRGP);
end
