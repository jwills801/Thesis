function plotElecRGPbySeaState()
% plotElecRGPbySeaState.m
% Line+marker plot of elecRGP by sea state for every drivetrain family:
% PassivePump/DHD2/DHD3/DHD4 (pressure-optimized, mAstar=5, from
% results/mAstarConvAndPressureOpt/pressureOptSummary.csv -- partial sea
% state coverage as of this writing) and EHA_mech/EHA_elec (full-length
% run, all 8 sea states, from results/EHA_allSeaStates/summary.csv).
% Missing (family, sea state) combinations are left as gaps (NaN), which
% simply break the line rather than leaving a misaligned/absent bar --
% unlike the grouped-bar version this replaces, this reads cleanly even
% though DHD3/DHD4 only have 4/8 and 1/8 sea states so far. Families with
% full 8/8 coverage are drawn as solid lines; partial-coverage families
% are dashed, to flag at a glance that they're still filling in (DHD4's
% single point renders as an isolated marker, no line). EHA_mech is
% de-emphasized (thin gray dashed line) since it's not a real contender
% for best electrical drivetrain -- it's here only to show how much
% mech-only optimization costs elecRGP relative to EHA_elec.
%
% Reflects the EHA copper-loss sign(T_Act) fix -- see diagnostics/ReadMe.md's
% "EHA copper-loss sign(T_Act) removed" section.
%
% Writes results/elecRGPbySeaState_withEHA_plot.png.
%
% Calls: none
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));

pressureOptFile = fullfile(repoRoot,'results','mAstarConvAndPressureOpt','pressureOptSummary.csv');
ehaFile = fullfile(repoRoot,'results','EHA_allSeaStates','summary.csv');

pressureOpt = readtable(pressureOptFile);
eha = readtable(ehaFile);

% pressureOptSummary.csv was generated while evaluate.m's DHD/PassivePump
% main-motor efficiency was temporarily 90% (since reverted back to 85%,
% see diagnostics/ReadMe.md's "DHD/PassivePump main-motor efficiency
% raised 85% -> 90%" section) -- rescale here rather than rerunning the
% whole pressure-optimization sweep, since this is an exact flat
% multiplicative correction (0.85/0.90) independent of sea state/family.
mainMotorRescale = 0.85/0.90;
pressureOpt.ElecRGP = pressureOpt.ElecRGP * mainMotorRescale;
pressureOpt.AveElecPow_kW = pressureOpt.AveElecPow_kW * mainMotorRescale;

nSeaStates = 8;
families = {'PassivePump','DHD2','DHD3','DHD4','EHA_mech','EHA_elec'};
nFamilies = numel(families);

elecRGP = nan(nSeaStates, nFamilies);
counts = zeros(1, nFamilies);

for f = 1:4 % PassivePump, DHD2, DHD3, DHD4
    rows = strcmp(pressureOpt.Family, families{f});
    idx = pressureOpt.SeaStateIdx(rows);
    elecRGP(idx, f) = pressureOpt.ElecRGP(rows);
    counts(f) = nnz(rows);
end
for f = 5:6 % EHA_mech, EHA_elec
    rows = strcmp(eha.Variant, families{f});
    idx = eha.SeaStateIdx(rows);
    elecRGP(idx, f) = eha.ElecRGP(rows);
    counts(f) = nnz(rows);
end

legendEntries = cell(1, nFamilies);
for f = 1:nFamilies
    legendEntries{f} = sprintf('%s (%d/%d pts)', families{f}, counts(f), nSeaStates);
end

% Distinct color per family (kept close to the old bar chart's palette);
% EHA_mech overridden to gray below so it recedes rather than competing
% for attention with the electrical-drivetrain comparison.
colors = [0.00 0.45 0.74;   % PassivePump - blue
          0.85 0.33 0.10;   % DHD2 - orange/red
          0.93 0.69 0.13;   % DHD3 - yellow
          0.49 0.18 0.56;   % DHD4 - purple
          0.50 0.50 0.50;   % EHA_mech - gray (de-emphasized)
          0.30 0.75 0.93];  % EHA_elec - cyan

fig = figure('Visible','off','Position',[1000 100 900 600]);
hold on
for f = 1:nFamilies
    fullCoverage = counts(f) == nSeaStates;
    if strcmp(families{f}, 'EHA_mech')
        lineWidth = 1.0; markerSize = 4; lineStyle = '--';
    else
        lineWidth = 2.0; markerSize = 6;
        if fullCoverage, lineStyle = '-'; else, lineStyle = '--'; end
    end
    plot(1:nSeaStates, elecRGP(:,f), 'o', 'LineStyle', lineStyle, ...
        'Color', colors(f,:), 'MarkerFaceColor', colors(f,:), ...
        'LineWidth', lineWidth, 'MarkerSize', markerSize);
end
hold off
xlabel('Sea state index')
ylabel('elecRGP')
xlim([0.5, nSeaStates+0.5])
xticks(1:nSeaStates)
title({'elecRGP by sea state, all families', ...
    '(DHD/PassivePump: pressure-optimized, mAstar=5; EHA: full-length run, post copper-loss-sign fix)', ...
    '(dashed = partial sea-state coverage so far; EHA_{mech} de-emphasized, not a real contender)'}, ...
    'FontSize', 11)
legend(legendEntries, 'Location', 'eastoutside')
grid on

resultsDir = fullfile(repoRoot,'results');
outFile = fullfile(resultsDir, 'elecRGPbySeaState_withEHA_plot.png');
set(fig,'renderer','painters'); % headless-safe rendering
exportgraphics(fig, outFile, 'Resolution', 150);
close(fig);
fprintf('Wrote %s\n', outFile);
end
