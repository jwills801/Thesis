function plotAstarHorizonResults()
% plotAstarHorizonResults.m
% Plots A* horizon length (m_Astar) vs. average electrical power, one
% line per sea state, from results/astarHorizon/summary.csv. Saves
% results/astarHorizon/plots/astarHorizonSweep.png.
% Calls: none
% Called by: Run_AstarHorizonSearch.m (automatically); fine standalone.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','astarHorizon');
plotDir = fullfile(resultsDir,'plots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

summaryFile = fullfile(resultsDir,'summary.csv');
if ~isfile(summaryFile)
    warning('plotAstarHorizonResults:noSummary','No summary.csv yet, nothing to plot.');
    return
end
summary = readtable(summaryFile);

seaStateIdxs = unique(summary.SeaStateIdx);
colors = lines(numel(seaStateIdxs));

fig = figure('Visible','off','Position',[100 100 700 500]);
hold on
legendEntries = {};
for k = 1:numel(seaStateIdxs)
    s = seaStateIdxs(k);
    rows = summary(summary.SeaStateIdx==s,:);
    rows = sortrows(rows,'mAstar');
    plot(rows.mAstar, rows.AveElecPow_kW, 'o-', 'Color', colors(k,:), ...
        'LineWidth', 1.5, 'MarkerFaceColor', colors(k,:));
    legendEntries{end+1} = sprintf('Sea state %d (Hs=%.2fm)', s, rows.Hs_m(1)); %#ok<AGROW>
end
xlabel('A* horizon length (m_{Astar})')
ylabel('Average Electrical Power [kW]')
title('DHD 2-rail: A* horizon length vs. electrical power')
legend(legendEntries, 'Location', 'best')
grid on
hold off

outFile = fullfile(plotDir,'astarHorizonSweep.png');
set(fig,'renderer','painters');
exportgraphics(fig, outFile, 'Resolution', 150);
close(fig);
fprintf('Wrote %s\n', outFile);
end
