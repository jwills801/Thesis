function plotAstarConvergenceResults()
% plotAstarConvergenceResults.m
% Plots results/astarConvergence100ms/summary.csv (from
% Run_AstarConvergenceStudy_100ms.m): elecRGP and mechRGP vs. m_Astar,
% one line per DHD rail count, for the single most-annual-energy sea
% state. Safe to run at any time, including mid-sweep.
%
% Calls: none
% Called by: none (run standalone once the study has some/all results)

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','astarConvergence100ms');
summaryFile = fullfile(resultsDir,'summary.csv');
if ~isfile(summaryFile)
    error('plotAstarConvergenceResults:noSummary','%s not found -- run Run_AstarConvergenceStudy_100ms.m first.', summaryFile);
end
T = readtable(summaryFile);

plotsDir = fullfile(resultsDir,'plots');
if ~exist(plotsDir,'dir'), mkdir(plotsDir); end

families = unique(T.Family);
colors = lines(numel(families));

fig = figure('Visible','off');
subplot(2,1,1); hold on
for f = 1:numel(families)
    rows = strcmp(T.Family, families{f});
    Tf = sortrows(T(rows,:),'MAstar');
    plot(Tf.MAstar*0.1, Tf.ElecRGP, '-o', 'Color', colors(f,:), 'DisplayName', families{f});
end
xlabel('Foresight time = m_{Astar} x 0.1s'); ylabel('Electrical RGP');
title('A* horizon convergence -- sea state 7 (most annual energy)');
legend('Location','best'); grid on

subplot(2,1,2); hold on
for f = 1:numel(families)
    rows = strcmp(T.Family, families{f});
    Tf = sortrows(T(rows,:),'MAstar');
    plot(Tf.MAstar*0.1, Tf.MechRGP, '-o', 'Color', colors(f,:), 'DisplayName', families{f});
end
xlabel('Foresight time = m_{Astar} x 0.1s'); ylabel('Mechanical RGP');
legend('Location','best'); grid on

set(fig,'renderer','painters');
outFile = fullfile(plotsDir,'astarConvergence.png');
saveas(fig, outFile);
close(fig);
fprintf('Wrote %s\n', outFile);
end
