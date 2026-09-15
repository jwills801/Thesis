% sizeCylinderArea.m
% Iteratively resizes a drivetrain family's cylinder area so the literal
% peak hydraulic pressure reached in a full closed-loop run at a given
% (sizing) sea state equals targetPressure, holding that pressure ceiling
% fixed throughout. For DHD, regenerates the switch-loss map fresh at
% each area iteration (cheap here since pressureRails doesn't change
% during sizing, unlike the later per-sea-state highPressure sweep).
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, makeSwitchLossMap.m
% Called by: none yet within this repo (invoked ad hoc from the
%   drivetrain-comparison sizing script)
function [capArea,rodArea,history] = sizeCylinderArea(runParamsTemplate,seaState,targetPressure,maxIter)
% Iteratively resizes the cylinder area for a drivetrain family so the
% LITERAL PEAK hydraulic pressure reached during a full closed-loop run
% at `seaState` equals `targetPressure` (Pa). highPressure (the ceiling
% rail, for PP/DHD) is held fixed at targetPressure throughout -- only
% area is being solved for. Ratio convention: EHA keeps capArea=rodArea
% (symmetric); PassivePump/DHD keep capArea=1.5*rodArea.
%
% For DHD, the switch-loss map only depends on hyd.pressureRails (a
% function of highPressure, fixed here) and capArea/rodArea (fixed at the
% iteration's current guess) -- NOT on the outcome we're converging on --
% so it's regenerated once per iteration's actual area, not per sea state
% grid point (that expensive per-grid-point regeneration only matters
% later when highPressure itself is swept, which optimizePressure.m
% already does correctly).
if nargin<4, maxIter=8; end

capArea = runParamsTemplate.capArea; % initial guess
history = struct('capArea',{},'peakPressure',{});

for iter=1:maxIter
    runParams = runParamsTemplate;
    if strcmp(runParams.drive,'EHA')
        runParams.capArea = capArea; runParams.rodArea = capArea;
    else
        runParams.capArea = capArea; runParams.rodArea = capArea/1.5;
    end

    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = seaState.Hs;
    params.simu.peakPeriod = seaState.Tp;

    if strcmp(runParams.drive,'DHD')
        params.hyd.switchMap = makeSwitchLossMap(params.hyd);
    end

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);

    F = dyn.u ./ params.hyd.Force2Torque(dyn.theta);
    deltaP = F/params.hyd.capArea;
    peakPressure = max(abs(deltaP));

    history(iter).capArea = capArea;
    history(iter).peakPressure = peakPressure;
    fprintf('    iter %d: capArea=%.5f m^2  ->  peakPressure=%.2f MPa\n', iter, capArea, peakPressure/1e6);

    ratio = peakPressure/targetPressure;
    if abs(ratio-1) < 0.02
        break
    end
    capArea = capArea*ratio;
end

if strcmp(runParamsTemplate.drive,'EHA')
    rodArea = capArea;
else
    rodArea = capArea/1.5;
end
end
