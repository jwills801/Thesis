function plotPressureSweeps()
% plotPressureSweeps.m
% For each sea state, plots the pressure grid search curve (candidate
% pressure vs. average electrical power) for every pressure-optimized
% drivetrain family (PassivePump, DHD2/3/4) on one figure, using the
% full grid data saved in each task's .fullGrid field (see
% Run_GridSearch.m). One PNG per sea state, saved to
% results/gridSearch/plots/ -- safe to run any time task_*.mat files with
% needsPressureOpt data exist, including mid-sweep.
%
% Calls: none
% Called by: Run_GridSearch.m (automatically, after aggregation)

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','gridSearch');
plotDir = fullfile(resultsDir,'plots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

files = dir(fullfile(resultsDir,'task_*.mat'));
pressureFamilies = {'PassivePump','DHD2','DHD3','DHD4'};
colors = lines(numel(pressureFamilies));

% Group tasks by sea state index
bySeaState = containers.Map('KeyType','double','ValueType','any');
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    if ~isfield(o,'fullGrid') || ~any(strcmp(o.label,pressureFamilies))
        continue % EHA tasks have no pressure grid to plot
    end
    if ~isKey(bySeaState,o.seaStateIdx)
        bySeaState(o.seaStateIdx) = struct('Hs',o.Hs,'Tp',o.Tp,'entries',{{}});
    end
    entry = bySeaState(o.seaStateIdx);
    entry.entries{end+1} = o;
    bySeaState(o.seaStateIdx) = entry;
end

seaStateIdxs = cell2mat(keys(bySeaState));
for s = sort(seaStateIdxs)
    entry = bySeaState(s);
    fig = figure('Visible','off','Position',[100 100 700 500]);
    hold on
    legendEntries = {};
    for e = 1:numel(entry.entries)
        o = entry.entries{e};
        famIdx = find(strcmp(o.label,pressureFamilies));
        plot(o.fullGrid.pressureGrid/1e6, o.fullGrid.aveElecPow/1e3, ...
            'o-', 'Color', colors(famIdx,:), 'LineWidth', 1.5, 'MarkerFaceColor', colors(famIdx,:));
        legendEntries{end+1} = o.label; %#ok<AGROW>
    end
    xlabel('Candidate highPressure [MPa]')
    ylabel('Average Electrical Power [kW]')
    title(sprintf('Sea state %d: Hs=%.2fm, Tp=%.2fs', s, entry.Hs, entry.Tp))
    legend(legendEntries, 'Location', 'best')
    grid on
    hold off

    outFile = fullfile(plotDir, sprintf('pressureSweep_seaState%02d.png', s));
    set(fig,'renderer','painters'); % headless-safe rendering
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
end
fprintf('Wrote %d pressure-sweep plot(s) to %s\n', numel(seaStateIdxs), plotDir);
end
