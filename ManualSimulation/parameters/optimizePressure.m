% optimizePressure.m
% Grid-searches runParams.highPressure over pressureGrid for a given sea
% state, running the full closed-loop simulation at each point, returning
% whichever maximizes electrical RGP. For DHD, regenerates (or reuses a
% passed-in denseSwitchMap for) the switch-loss map per point, since it's
% only valid for the exact rail pressures it was built at.
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m,
%   makeSwitchLossMap.m (unless denseSwitchMap is given)
% Called by: diagnostics/validatePhase2Subset.m
function result = optimizePressure(runParams, seaState, pressureGrid, simOverrides, denseSwitchMap)
% simOverrides (optional): struct of fields to overwrite on params.simu
% after getParameters (e.g. finalTime/rampTime/time, for a fast/short
% validation sweep). Omit for a full-length production run.
%
% denseSwitchMap (optional, DHD only): a switch-loss map already built
% over a dense pressure grid spanning the full range this sweep's
% candidate highPressure values fall within (see buildDenseSwitchMap.m).
% When given, it's reused as-is at every grid point instead of
% regenerating a fresh exact-rails map per point -- valid because
% consumers (getValveLoss.m etc.) already interpolate via interpn against
% whatever switchMap.PR axis is given, not against the actual rails in
% use. Omit to keep the original per-point regeneration behavior.
if nargin < 4
    simOverrides = struct();
end
if nargin < 5
    denseSwitchMap = [];
end
% Grid-searches runParams.highPressure (Pa) over pressureGrid for a given
% sea state (seaState.Hs [m], seaState.Tp [s]), running the full
% drivetrain simulation at each grid point and returning the pressure
% that maximizes electrical RGP, plus the full sweep curve (so
% unimodality can be checked visually rather than assumed -- this
% replaces the disabled secant-search OptPressure() in Run_All_Cases.m).
%
% Applies to rail-scheduled drivetrains: PassivePump/CoulombDamping and
% DHD (any controller with pressure rails, e.g. MPC_Astar). EHA has no
% pressure rails and should not be passed here.
%
% IMPORTANT -- DHD switching-loss map: parameters/getHydraulic.m normally
% loads a single precomputed parameters/SwitchMap.mat regardless of the
% current highPressure. That map is only physically valid for the exact
% rail pressures it was built at (makeSwitchLossMap.m bakes hyd.pressureRails
% into its PR grid). Sweeping highPressure away from that value means the
% switching-loss interpolation in getValveLoss.m/getSwitchingLoss.m gets
% evaluated off its native grid -- silently wrong losses, not an error.
% This function regenerates the switch-loss map fresh at each grid point
% for DHD so the sweep is physically consistent. That regeneration
% (makeSwitchLossMap.m) is the expensive part of this function -- see the
% ReadMe in diagnostics/ for measured cost per call before running a large
% sweep.

result = struct();
result.pressureGrid = pressureGrid(:);
result.mechRGP = NaN(size(result.pressureGrid));
result.elecRGP = NaN(size(result.pressureGrid));

for i = 1:numel(pressureGrid)
    thisRunParams = runParams;
    thisRunParams.highPressure = pressureGrid(i);

    params = getParameters(thisRunParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = seaState.Hs;
    params.simu.peakPeriod = seaState.Tp;
    overrideFields = fieldnames(simOverrides);
    for f = 1:numel(overrideFields)
        params.simu.(overrideFields{f}) = simOverrides.(overrideFields{f});
    end

    if strcmp(runParams.drive,'DHD')
        if ~isempty(denseSwitchMap)
            params.hyd.switchMap = denseSwitchMap;
        else
            params.hyd.switchMap = makeSwitchLossMap(params.hyd);
        end
    end

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    result.mechRGP(i) = ev.mechRGP;
    result.elecRGP(i) = ev.elecRGP;
end

[~,bestInd] = max(result.elecRGP);
result.bestPressure = result.pressureGrid(bestInd);
result.bestElecRGP = result.elecRGP(bestInd);
end
