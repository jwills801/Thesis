function plotEHAFixedVsVariableBySeaState()
% plotEHAFixedVsVariableBySeaState.m
% Line+marker plot (same style as plotElecRGPbySeaState.m) of mechRGP and
% elecRGP by sea state for the original large-displacement (Scale x1,
% real 107cc/rev pump) EHA fixed- vs. variable-displacement comparison --
% results/EHA_fixedVsVariable/summary.csv, from Run_EHA_FixedVsVariable.m.
% All 4/4 EHA variants have full 8/8 sea-state coverage, so no dashed/
% partial-coverage handling is needed here (unlike the DHD/PassivePump
% plot). *_mech (mechanical-only objective) is de-emphasized (thin gray)
% since it's not a real contender, same convention as
% plotElecRGPbySeaState.m -- included only to show what loss-aware
% control buys over mechanical-only.
%
% Writes results/EHA_fixedVsVariable_bySeaState_plot.png.
%
% Calls: none
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
T = readtable(fullfile(repoRoot,'results','EHA_fixedVsVariable','summary.csv'));

nSeaStates = 8;
variants = {'EHA_fixed_elec','EHA_var_elec','EHA_fixed_mech','EHA_var_mech'};
nVariants = numel(variants);

mechRGP = nan(nSeaStates,nVariants);
elecRGP = nan(nSeaStates,nVariants);
for v = 1:nVariants
    rows = strcmp(T.Variant, variants{v});
    idx = T.SeaStateIdx(rows);
    mechRGP(idx,v) = T.MechRGP(rows);
    elecRGP(idx,v) = T.ElecRGP(rows);
end

% color, style, width, marker size per variant
styles = struct( ...
    'EHA_fixed_elec', struct('color',[0.00 0.45 0.74],'style','-','width',2.0,'msize',6), ...
    'EHA_var_elec',   struct('color',[0.85 0.33 0.10],'style','-','width',2.0,'msize',6), ...
    'EHA_fixed_mech', struct('color',[0.30 0.60 0.90],'style','--','width',1.0,'msize',4), ...
    'EHA_var_mech',   struct('color',[0.93 0.60 0.45],'style','--','width',1.0,'msize',4) ...
    );

fig = figure('Visible','off','Position',[100 100 900 750]);

subplot(2,1,1)
hold on
for v = 1:nVariants
    s = styles.(variants{v});
    plot(1:nSeaStates, mechRGP(:,v), 'o', 'LineStyle', s.style, 'Color', s.color, ...
        'MarkerFaceColor', s.color, 'LineWidth', s.width, 'MarkerSize', s.msize);
end
hold off
ylabel('mechRGP')
xlim([0.5, nSeaStates+0.5]); xticks(1:nSeaStates)
title('EHA fixed- vs. variable-displacement, Scale x1 (real 107cc/rev pump)')
grid on
legend(strrep(variants,'_','\_'), 'Location', 'eastoutside')

subplot(2,1,2)
hold on
for v = 1:nVariants
    s = styles.(variants{v});
    plot(1:nSeaStates, elecRGP(:,v), 'o', 'LineStyle', s.style, 'Color', s.color, ...
        'MarkerFaceColor', s.color, 'LineWidth', s.width, 'MarkerSize', s.msize);
end
yline(0,'k-','LineWidth',0.75);
hold off
xlabel('Sea state index'); ylabel('elecRGP')
xlim([0.5, nSeaStates+0.5]); xticks(1:nSeaStates)
grid on
legend(strrep(variants,'_','\_'), 'Location', 'eastoutside')

resultsDir = fullfile(repoRoot,'results');
outFile = fullfile(resultsDir, 'EHA_fixedVsVariable_bySeaState_plot.png');
set(fig,'renderer','painters');
exportgraphics(fig, outFile, 'Resolution', 150);
close(fig);
fprintf('Wrote %s\n', outFile);
end
