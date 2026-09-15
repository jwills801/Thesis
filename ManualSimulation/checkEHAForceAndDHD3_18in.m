% checkEHAForceAndDHD3_18in.m
% One-off combined script:
%  1) Pulls the peak-force trajectory from an EHA sizing-equivalent run
%     (same sizing sea state, sized capArea) to check whether EHA's ~2.5x
%     bigger peak force (vs PassivePump/DHD) is a brief spike or sustained.
%  2) Builds a fresh 100ms dense switch-loss map for an 18in cap bore
%     (capArea=0.164173 m^2, rodArea=0.109449 m^2 -- exactly 1.5x the
%     current 14.70in DHD family's areas; note the new rodArea exactly
%     equals the CURRENT capArea).
%  3) Runs DHD3 at sea state 8 (Hs=4.995m, Tp=12.725m -- the sizing sea
%     state, confirmed as "most energetic" for this test) at both the
%     current 14.70in bore and the new 18in bore, same controlDT=0.1,
%     mAstar=12, astarIterMax=100000, highPressure=35e6, for a direct
%     side-by-side comparison.
%
% Requires results/denseSwitchMaps_100ms.mat already rebuilt with the
% corrected valve dynamics (rebuildSwitchMaps.m) before this runs.
%
% Calls: optimization/buildDenseSwitchMap.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'), fullfile(repoRoot,'wave'), ...
    fullfile(repoRoot,'control'), fullfile(repoRoot,'dynamics'), ...
    fullfile(repoRoot,'evaluation'));

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));

Hs8 = 4.99525316455696;
Tp8 = 12.7246835443038;

%% (1) EHA peak-force trajectory check
fprintf('=== 1) EHA peak-force trajectory at the sizing sea state ===\n');
runParamsEHA = struct('drive','EHA','controller','MPC_QP','considerLosses',0, ...
    'capArea',S.sizedAreas.EHA.capArea,'rodArea',S.sizedAreas.EHA.rodArea, ...
    'ehaFixedDisplacement',true,'shaftInertia',249.9);
paramsEHA = getParameters(runParamsEHA);
paramsEHA.simu.makePlots = false;
paramsEHA.simu.sigWaveHeight = Hs8;
paramsEHA.simu.peakPeriod = Tp8;
waveEHA = generateExcitingTorque(paramsEHA);
ctrlEHA = getControl(paramsEHA,waveEHA);
dynEHA = timeLoop(paramsEHA,waveEHA,ctrlEHA);

F = dynEHA.u ./ paramsEHA.hyd.Force2Torque(dynEHA.theta);
peakPressure = max(abs(F))/paramsEHA.hyd.capArea;
fprintf('Peak force: %.3f MN, implied peak pressure: %.2f MPa (sizing target 35MPa)\n', max(abs(F))/1e6, peakPressure/1e6);

rampInd = round(paramsEHA.simu.rampTime/paramsEHA.simu.dt);
Fpost = F(rampInd:end);
threshold90 = 0.9*max(abs(Fpost));
fracAbove90 = mean(abs(Fpost) > threshold90);
fracAbove50 = mean(abs(Fpost) > 0.5*max(abs(Fpost)));
fprintf('Fraction of post-ramp time with |F| > 90%% of peak: %.3f%%\n', 100*fracAbove90);
fprintf('Fraction of post-ramp time with |F| > 50%% of peak: %.3f%%\n', 100*fracAbove50);
fprintf('RMS force (post-ramp): %.3f MN (peak/RMS ratio: %.2f)\n', rms(Fpost)/1e6, max(abs(Fpost))/rms(Fpost));

%% (2) Build 18in-bore 100ms switch map
fprintf('\n=== 2) Building 18in-bore 100ms switch map ===\n');
capArea18 = 0.164173; rodArea18 = 0.109449;
hydTmp18 = struct('capArea',capArea18,'rodArea',rodArea18,'stroke',5);
tic;
denseMap18 = buildDenseSwitchMap(hydTmp18, 0.5e6, 35e6, 10, 0.1);
fprintf('Build took %.1f minutes\n', toc/60);
save(fullfile(repoRoot,'results','denseSwitchMap_18in_100ms.mat'),'denseMap18');

%% (3) DHD3 comparison: current bore vs 18in bore, sea state 8, 100ms, mAstar=12
D100 = load(fullfile(repoRoot,'results','denseSwitchMaps_100ms.mat'));
denseMapCurrent = D100.denseMapDHD2; % valid for all rail counts, see its header comment

fprintf('\n=== 3) DHD3 @ sea state 8, 100ms, mAstar=12: current (14.70in) vs 18in bore ===\n');

configs = {
    'current (14.70in)', S.sizedAreas.DHD3.capArea, S.sizedAreas.DHD3.rodArea, denseMapCurrent
    '18in', capArea18, rodArea18, denseMap18
};

for i = 1:size(configs,1)
    label = configs{i,1}; cA = configs{i,2}; rA = configs{i,3}; map = configs{i,4};
    runParams = struct('drive','DHD','controller','MPC_Astar','pressure_rails',3, ...
        'considerLosses',1,'capArea',cA,'rodArea',rA, ...
        'highPressure',35e6,'astarIterMax',100000,'controlDT',0.1,'mAstar',12);
    params = getParameters(runParams);
    params.simu.makePlots = false;
    params.simu.sigWaveHeight = Hs8;
    params.simu.peakPeriod = Tp8;
    params.hyd.switchMap = map;

    wave = generateExcitingTorque(params);
    ctrl = getControl(params,wave);
    dyn = timeLoop(params,wave,ctrl);
    ev = evaluate(params,dyn,ctrl);

    fprintf('%-20s capArea=%.6f mechRGP=%.4f elecRGP=%.4f aveMechPow=%.1fkW aveElecPow=%.1fkW capHits=%d\n', ...
        label, cA, ev.mechRGP, ev.elecRGP, ev.aveMechPow/1e3, ev.aveElecPow/1e3, ev.nAstarCapHits);
end
