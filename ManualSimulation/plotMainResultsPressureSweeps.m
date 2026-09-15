function plotMainResultsPressureSweeps()
% plotMainResultsPressureSweeps.m
% For each rail-scheduled family (PassivePump, DHD2, DHD3, DHD4), one
% figure with an 8-panel (2x4) grid, one panel per sea state, showing
% elecRGP vs. pressure for every task_<family>_ss<N>_p*.mat found in
% results/mainResults/ -- i.e. the actual pressure-sweep curve, not just
% the single selected best point. The chosen point (matching
% results/mainResults.csv's selection logic -- position-constrained
% argmax for PassivePump, plain argmax for DHD) is marked with a gold
% star. Data-completeness and PassivePump's position-target status are
% intentionally NOT shown visually (per feedback); panels with no data
% yet are left blank with a note.
%
% Color/marker/line-style per family come from plotting/drivetrainStyle.m
% -- the shared convention used across every drivetrain-comparison plot
% in this repo.
%
% Writes results/mainResults_pressureSweep_<family>.png (one per family).
%
% Calls: plotting/drivetrainStyle.m
% Called by: none (top-level entry point, rerun anytime against whatever
%   task files exist on disk)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'plotting'));
resultsDir = fullfile(repoRoot,'results','mainResults');

nSeaStates = 8;
POSITION_TARGET_DEG = 5;
families = {'PassivePump','DHD2','DHD3','DHD4'};

for f = 1:numel(families)
    family = families{f};
    s = drivetrainStyle(family);
    fig = figure('Visible','off','Position',[100 100 1400 700]);

    for ss = 1:nSeaStates
        subplot(2,4,ss)
        pattern = fullfile(resultsDir, sprintf('task_%s_ss%d_p*.mat', family, ss));
        files = dir(pattern);

        if isempty(files)
            text(0.5,0.5,'no data yet','HorizontalAlignment','center');
            axis off
            title(sprintf('SS%d', ss));
            continue
        end

        pressures = nan(numel(files),1); elecRGP = nan(numel(files),1);
        maxThetas = nan(numel(files),1); minThetas = nan(numel(files),1);
        for k = 1:numel(files)
            L = load(fullfile(files(k).folder, files(k).name));
            pressures(k) = L.out.pressure/1e6;
            elecRGP(k) = L.out.elecRGP;
            maxThetas(k) = L.out.maxThetaDeg;
            minThetas(k) = L.out.minThetaDeg;
        end
        [pressures, order] = sort(pressures); elecRGP = elecRGP(order);
        maxThetas = maxThetas(order); minThetas = minThetas(order);

        hold on
        plot(pressures, elecRGP, s.marker, 'LineStyle', s.lineStyle, 'Color', s.color, ...
            'MarkerFaceColor', s.color, 'MarkerEdgeColor', s.color, 'LineWidth', 1.5, 'MarkerSize', 7);

        % Symmetric +/-5deg position target, applied to every family:
        % among points clearing it on both sides, pick argmax(elecRGP);
        % otherwise fall back to whichever already-run pressure is
        % closest to 5MPa -- never a new simulation.
        candidates = find(maxThetas > POSITION_TARGET_DEG & minThetas < -POSITION_TARGET_DEG);
        if isempty(candidates)
            [~,bestInd] = min(abs(pressures - 5));
        else
            [~,bestLocal] = max(elecRGP(candidates)); bestInd = candidates(bestLocal);
        end
        plot(pressures(bestInd), elecRGP(bestInd), 'p', 'MarkerFaceColor', [0.93 0.69 0.13], ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 14);
        hold off
        grid on
        xlabel('Pressure (MPa)'); ylabel('elecRGP');
        title(sprintf('SS%d (%d pts)', ss, numel(files)));
    end

    sgtitle(sprintf('%s pressure sweep by sea state (star = selected)', family));

    outFile = fullfile(repoRoot,'results',sprintf('mainResults_pressureSweep_%s.png', family));
    set(fig,'renderer','painters');
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
    fprintf('Wrote %s\n', outFile);
end
end
