% Run_EHASpeedCheck.m
% Checks whether the EHA sizing pipeline's assumed design-point piston
% velocity (vMax=1 m/s, parameters/getHydraulic.m:82 -- used to derive
% maxFlow=vMax*capArea, which in turn sets Scale, the EHA pump's total
% installed displacement, see models/makeEHALossMap*.m) is actually
% reached in real closed-loop runs, across all 8 Humboldt sea states, for
% both EHA_fixed_elec and EHA_var_elec (Scale x1, real 107cc/rev pump,
% same sizing as results/EHA_fixedVsVariable/summary.csv). Also reports
% the flap's actual angular speed (thetaDot) and position (theta)
% magnitudes reached, since those drive the piston velocity via
% hyd.dLdt(theta,thetaDot).
%
% Calls: parameters/getParameters.m, wave/generateExcitingTorque.m,
%   control/getControl.m, dynamics/timeLoop.m, evaluation/evaluate.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'wave'), fullfile(repoRoot,'control'), ...
    fullfile(repoRoot,'dynamics'), fullfile(repoRoot,'evaluation'));

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
rodArea = S.sizedAreas.EHA.rodArea;

seaStates = humboldtSeaStates();
VMAX_DESIGN = 1; % m/s -- parameters/getHydraulic.m's sizing assumption

SEA_STATES_TO_RUN = 8; % SS8 only -- the sizing sea state (35MPa peak pressure target)
cases = {'EHA_var_elec', false; 'EHA_fixed_elec', true}; % elec-aware (considerLosses=1) only -- no mech-only EHA cases

fprintf('%-16s %2s  %8s  %10s  %10s  %10s  %10s\n', 'Variant','SS','Hs(m)','peakVel(m/s)','peakThDot','maxTheta(deg)','minTheta(deg)');
results = struct('label',{},'ss',{},'peakVel',{},'peakThetaDot',{},'maxThetaDeg',{},'minThetaDeg',{});
for c = 1:size(cases,1)
    label = cases{c,1}; fixedDisp = cases{c,2};
    for ss = SEA_STATES_TO_RUN
        Hs = seaStates.Hs(ss); Tp = seaStates.Tp(ss);

        runParams = struct('drive','EHA','controller','MPC_QP','considerLosses',1, ...
            'capArea',capArea,'rodArea',rodArea,'controlDT',0.1,'ehaFixedDisplacement',fixedDisp);
        if fixedDisp
            runParams.shaftInertia = 249.9;
        end

        params = getParameters(runParams);
        params.simu.makePlots = false;
        params.simu.sigWaveHeight = Hs; params.simu.peakPeriod = Tp;

        wave = generateExcitingTorque(params);
        ctrl = getControl(params,wave);
        dyn = timeLoop(params,wave,ctrl);

        thetaDot = dyn.thetaDot; theta = dyn.theta;
        vel = params.hyd.dLdt(theta,thetaDot); % actual piston velocity, m/s

        rampInd = round(params.simu.rampTime/params.simu.dt);
        postRamp = rampInd:length(dyn.t);

        peakVel = max(abs(vel(postRamp)));
        peakThetaDot = max(abs(thetaDot(postRamp)));
        maxThetaDeg = max(theta(postRamp))*180/pi;
        minThetaDeg = min(theta(postRamp))*180/pi;

        fprintf('%-16s %2d  %8.2f  %10.3f  %10.4f  %10.2f  %10.2f\n', ...
            label, ss, Hs, peakVel, peakThetaDot, maxThetaDeg, minThetaDeg);

        results(end+1) = struct('label',label,'ss',ss,'peakVel',peakVel, ...
            'peakThetaDot',peakThetaDot,'maxThetaDeg',maxThetaDeg,'minThetaDeg',minThetaDeg); %#ok<AGROW>
    end
end

allPeakVel = [results.peakVel];
[worstVel,worstInd] = max(allPeakVel);
fprintf('\nDesign assumption: vMax = %.1f m/s\n', VMAX_DESIGN);
fprintf('Actual peak piston velocity across all runs: %.3f m/s (%s, sea state %d) -- %.0f%% of design vMax\n', ...
    worstVel, results(worstInd).label, results(worstInd).ss, 100*worstVel/VMAX_DESIGN);

T = struct2table(results);
writetable(T, fullfile(repoRoot,'results','EHA_speedCheck.csv'));
fprintf('Wrote results/EHA_speedCheck.csv\n');
