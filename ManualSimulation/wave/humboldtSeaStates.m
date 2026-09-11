function seaStates = humboldtSeaStates(makePlot)
% Reads humboldtManualBinsCells.csv and collapses each of the 8
% manually-assigned bins to a single occurrence-weighted representative
% (Hs,Tp) sea state, for the drivetrain comparison sweep. Uses Te (energy
% period) directly as Tp -- flagged, not converted (Tp=Te/0.858 would be
% the standard PM conversion).
%
% Input CSV columns: Bin, Hm0_min_m, Hm0_max_m, Te_min_s, Te_max_s,
% Percent_Occurrence, Percent_Total_Energy. Bin membership is fixed by
% hand (not re-clustered).
%
% makePlot (optional, default false): if true, also plots the Hm0/Te JPD
% heatmap (occurrence and energy panels) with bin boundaries and
% representative points overlaid, and writes humboldtBinSummary.csv.
% Off by default so routine callers (e.g. diagnostics/validatePhase2Subset.m)
% don't pop a figure every call.
%
% Calls: none
% Called by: diagnostics/validatePhase2Subset.m
if nargin < 1
    makePlot = false;
end

here = fileparts(mfilename('fullpath'));
csvPath = fullfile(here,'humboldtManualBinsCells.csv');
T = readtable(csvPath);

Hm0c = (T.Hm0_min_m + T.Hm0_max_m)/2;
Tec  = (T.Te_min_s  + T.Te_max_s)/2;
occ  = T.Percent_Occurrence;
eng  = T.Percent_Total_Energy;
occ(isnan(occ)) = 0;
eng(isnan(eng)) = 0;
binId = T.Bin;
nBins = max(binId);

Hs = NaN(nBins,1); Tp = NaN(nBins,1); probability = NaN(nBins,1);
for b = 1:nBins
    idx = binId == b;
    probability(b) = sum(occ(idx));
    Hs(b) = sum(Hm0c(idx).*occ(idx)) / sum(occ(idx));
    Tp(b) = sum(Tec(idx).*occ(idx)) / sum(occ(idx));
end

% Raw occurrence percentages sum to ~100 but not exactly; normalize to 1.
probability = probability/sum(probability);

seaStates = table(Hs, Tp, probability);

if ~makePlot
    return
end

%% ---------------- Per-bin summary table (for the plot + CSV export) ----------------
binSummary = table((1:nBins)', 'VariableNames', {'Bin'});
binSummary.N_cells        = zeros(nBins,1);
binSummary.Occurrence_pct = zeros(nBins,1);
binSummary.Energy_pct     = zeros(nBins,1);
binSummary.Hm0_rep_m      = zeros(nBins,1);
binSummary.Te_rep_s       = zeros(nBins,1);
binSummary.Hm0_range_m    = strings(nBins,1);
binSummary.Te_range_s     = strings(nBins,1);

for b = 1:nBins
    idx = binId == b;
    binSummary.N_cells(b)        = sum(idx);
    binSummary.Occurrence_pct(b) = sum(occ(idx));
    binSummary.Energy_pct(b)     = sum(eng(idx));
    binSummary.Hm0_rep_m(b) = Hs(b);
    binSummary.Te_rep_s(b)  = Tp(b);
    binSummary.Hm0_range_m(b) = sprintf('%.1f-%.1f', min(T.Hm0_min_m(idx)), max(T.Hm0_max_m(idx)));
    binSummary.Te_range_s(b)  = sprintf('%.0f-%.0f',  min(T.Te_min_s(idx)),  max(T.Te_max_s(idx)));
end

disp(binSummary)
writetable(binSummary, fullfile(here,'humboldtBinSummary.csv'));

%% ---------------- Build grids for plotting ----------------
hm0Edges = 0:0.5:9;
teEdges  = 3:1:18;
nH = length(hm0Edges) - 1;
nT = length(teEdges) - 1;

occGrid = nan(nH, nT);
engGrid = nan(nH, nT);
binGrid = nan(nH, nT);

for i = 1:height(T)
    hi = find(hm0Edges == T.Hm0_min_m(i), 1);
    ti = find(teEdges  == T.Te_min_s(i), 1);
    if isempty(hi) || isempty(ti), continue; end
    occGrid(hi, ti) = occ(i);
    engGrid(hi, ti) = eng(i);
    binGrid(hi, ti) = binId(i);
end

%% ---------------- Colormap (YlOrRd-style, punchy for heatmap emphasis) ----------------
ylOrRd = [255 255 204; 255 237 160; 254 217 118; 254 178 76; ...
          253 141 60;  252 78 42;   227 26 28;   189 0 38; 128 0 38] / 255;
cmap = interp1(1:size(ylOrRd,1), ylOrRd, linspace(1, size(ylOrRd,1), 256));

%% ---------------- Plot ----------------
figure('Position', [100 100 950 950], 'Color', 'w');

panelData   = {occGrid, engGrid};
panelTitles = ["Humboldt: % Occurrence", "Humboldt: % of Total Energy"];
cbLabels    = ["% Occurrence", "% of Total Energy"];
panelMax    = [6, 4.5];

for p = 1:2
    ax = subplot(2,1,p);
    gridData = panelData{p};

    Cpad = [gridData, nan(nH,1); nan(1, nT+1)];
    h = pcolor(teEdges, hm0Edges, Cpad);
    shading flat
    set(h, 'EdgeColor', [0.92 0.92 0.92], 'LineWidth', 0.3);
    colormap(ax, cmap);
    caxis(ax, [0 panelMax(p)]);
    cb = colorbar;
    ylabel(cb, cbLabels(p), 'FontSize', 9);
    hold on

    for i = 1:nH
        for j = 1:nT
            v = gridData(i,j);
            if ~isnan(v)
                frac = min(v / panelMax(p), 1);
                txtColor = 'k';
                if frac > 0.55, txtColor = 'w'; end
                text(teEdges(j)+0.5, hm0Edges(i)+0.25, sprintf('%.2f', v), ...
                    'HorizontalAlignment', 'center', 'FontSize', 6, 'Color', txtColor);
            end
        end
    end

    edgeColor = [0.25 0.25 0.25];
    for i = 1:nH
        for j = 1:nT
            b = binGrid(i,j);
            if isnan(b), continue; end
            if j == nT || isnan(binGrid(i,j+1)) || binGrid(i,j+1) ~= b
                plot([teEdges(j+1) teEdges(j+1)], [hm0Edges(i) hm0Edges(i+1)], 'Color', edgeColor, 'LineWidth', 0.9);
            end
            if j == 1 || isnan(binGrid(i,j-1)) || binGrid(i,j-1) ~= b
                plot([teEdges(j) teEdges(j)], [hm0Edges(i) hm0Edges(i+1)], 'Color', edgeColor, 'LineWidth', 0.9);
            end
            if i == nH || isnan(binGrid(i+1,j)) || binGrid(i+1,j) ~= b
                plot([teEdges(j) teEdges(j+1)], [hm0Edges(i+1) hm0Edges(i+1)], 'Color', edgeColor, 'LineWidth', 0.9);
            end
            if i == 1 || isnan(binGrid(i-1,j)) || binGrid(i-1,j) ~= b
                plot([teEdges(j) teEdges(j+1)], [hm0Edges(i) hm0Edges(i)], 'Color', edgeColor, 'LineWidth', 0.9);
            end
        end
    end

    for b = 1:nBins
        hrep = binSummary.Hm0_rep_m(b);
        trep = binSummary.Te_rep_s(b);
        plot(trep, hrep, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
        text(trep + 0.15, hrep + 0.15, num2str(b), 'FontSize', 8, ...
            'FontWeight', 'bold', 'Color', 'k', 'BackgroundColor', [1 1 1]);
    end

    xlabel('Energy Period, T_e [s]');
    ylabel('Significant Wave Height, H_{m0} [m]');
    title(panelTitles(p) + " | manually assigned bins");
    xlim([3 18]); ylim([0 9]);
    set(ax, 'Layer', 'top');
    hold off
end

sgtitle('Humboldt Bay JPD with Manually Assigned Sea-State Bins', ...
    'FontSize', 12, 'FontWeight', 'bold');

exportgraphics(gcf, fullfile(here,'humboldtJpdManualBins.png'), 'Resolution', 200);
end
