function plotMainResultsBySeaState()
% plotMainResultsBySeaState.m
% Line+marker plot (mechRGP top, elecRGP bottom) by sea state for every
% drivetrain family in results/mainResults.csv (PassivePump, DHD2, DHD3,
% DHD4, EHA_fixed_elec), each using its selected/best pressure (or single
% run, for EHA) per sea state. Missing sea states (no data yet) simply
% break the line -- same convention as plotElecRGPbySeaState.m. Data-
% completeness is intentionally NOT shown (per feedback -- not of
% interest here).
%
% Color/marker/line-style per family come from plotting/drivetrainStyle.m
% -- the shared convention used across every drivetrain-comparison plot
% in this repo, so this figure matches the pressure-sweep plots and any
% future ones.
%
% Writes results/mainResults_bySeaState_plot.png.
%
% Calls: plotting/drivetrainStyle.m
% Called by: none (top-level entry point, rerun anytime)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'plotting'));
T = readtable(fullfile(repoRoot,'results','mainResults.csv'));

nSeaStates = 8;
families = {'PassivePump','DHD2','DHD3','DHD4','EHA_fixed_elec'};
styleNames = {'PassivePump','DHD2','DHD3','DHD4','EHA'};

fig = figure('Visible','off','Position',[100 100 950 800]);

for panel = 1:2
    subplot(2,1,panel)
    hold on
    legendEntries = {};
    for f = 1:numel(families)
        family = families{f};
        fam = T(strcmp(T.Family,family),:);
        if isempty(fam), continue; end

        y = nan(nSeaStates,1);
        for i = 1:height(fam)
            ss = fam.SeaStateIdx(i);
            if panel == 1
                y(ss) = fam.MechRGP(i);
            else
                y(ss) = fam.ElecRGP(i);
            end
        end

        s = drivetrainStyle(styleNames{f});
        plot(1:nSeaStates, y, s.marker, 'LineStyle', s.lineStyle, 'Color', s.color, ...
            'MarkerFaceColor', s.color, 'MarkerEdgeColor', s.color, 'LineWidth', 2, 'MarkerSize', 7);
        legendEntries{end+1} = strrep(family,'_','\_'); %#ok<AGROW>
    end
    yline(0,'k-','LineWidth',0.75);
    hold off
    grid on
    xlim([0.5, nSeaStates+0.5]); xticks(1:nSeaStates);
    if panel == 1
        ylabel('mechRGP')
        title('Main results by sea state')
    else
        ylabel('elecRGP'); xlabel('Sea state index')
    end
    legend(legendEntries, 'Location', 'eastoutside')
end

outFile = fullfile(repoRoot,'results','mainResults_bySeaState_plot.png');
set(fig,'renderer','painters');
exportgraphics(fig, outFile, 'Resolution', 150);
close(fig);
fprintf('Wrote %s\n', outFile);
end
