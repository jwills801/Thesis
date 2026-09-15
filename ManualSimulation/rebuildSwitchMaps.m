% rebuildSwitchMaps.m
% Rebuilds both dense switch-loss maps (results/denseSwitchMaps.mat at
% 200ms, results/denseSwitchMaps_100ms.mat at 100ms) after the valve
% natural-frequency correction in models/makeSwitchLossMap.m (25Hz ->
% 20.11Hz, to actually match the valve's 26ms datasheet response time --
% the old value gave ~20.9ms). Every switch-loss-dependent DHD result
% computed before this rebuild used the wrong valve dynamics.
%
% ONE map per switchTime is built and reused for DHD2/DHD3/DHD4 --
% confirmed identical capArea/rodArea across rail counts (see
% optimization/buildDenseSwitchMap.m's comment) -- results/denseSwitchMaps.mat
% keeps its historical per-family struct shape (DHD2/DHD3/DHD4 fields, all
% pointing at the same map) only so Run_GridSearch.m doesn't need editing;
% results/denseSwitchMaps_100ms.mat keeps its existing single-map shape
% (denseMapDHD2, used generically -- see Run_AstarConvergenceStudy_100ms.m).
%
% Calls: optimization/buildDenseSwitchMap.m
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'), ...
    fullfile(repoRoot,'optimization'));

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
hydTmp = struct('capArea', S.sizedAreas.DHD2.capArea, 'rodArea', S.sizedAreas.DHD2.rodArea, 'stroke', 5);

fprintf('=== Building 200ms dense switch map (corrected valve dynamics) ===\n');
tic;
map200 = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.2);
fprintf('200ms build took %.1f minutes\n', toc/60);
denseSwitchMaps = struct('DHD2',map200,'DHD3',map200,'DHD4',map200);
save(fullfile(repoRoot,'results','denseSwitchMaps.mat'),'denseSwitchMaps');

fprintf('\n=== Building 100ms dense switch map (corrected valve dynamics) ===\n');
tic;
denseMapDHD2 = buildDenseSwitchMap(hydTmp, 0.5e6, 35e6, 10, 0.1);
fprintf('100ms build took %.1f minutes\n', toc/60);
save(fullfile(repoRoot,'results','denseSwitchMaps_100ms.mat'),'denseMapDHD2');

fprintf('\nBoth switch maps rebuilt with the corrected 20.11Hz valve dynamics.\n');
