function plotFinalPressureSweeps()
% plotFinalPressureSweeps.m
% For each of the 8 sea states, plots pressure vs. average electrical
% power for the 4 pressure-optimized families (PassivePump, DHD2/3/4)
% from results/gridSearchFinal/task_*.mat's saved fullGrid curves.
% Safe to run at any time once at least some tasks have completed.
%
% Calls: none
% Called by: Run_FinalGridSearch_100ms.m (automatically); also fine to
%   call directly/standalone at any time.

repoRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot,'results','gridSearchFinal');
plotsDir = fullfile(resultsDir,'plots');
if ~exist(plotsDir,'dir'), mkdir(plotsDir); end

files = dir(fullfile(resultsDir,'task_*.mat'));
if isempty(files)
    warning('plotFinalPressureSweeps:noTasks','No task_*.mat files found yet.');
    return
end

families = {'PassivePump','DHD2','DHD3','DHD4'};
colors = lines(numel(families));

data = struct();
for i = 1:numel(files)
    d = load(fullfile(files(i).folder,files(i).name));
    o = d.out;
    if ~isfield(o,'fullGrid'), continue; end % EHA tasks have no pressure sweep
    key = sprintf('ss%d', o.seaStateIdx);
    if ~isfield(data,key), data.(key) = struct('Hs',o.Hs,'Tp',o.Tp,'families',struct()); end
    data.(key).families.(o.label) = o.fullGrid;
end

seaStateKeys = fieldnames(data);
for k = 1:numel(seaStateKeys)
    key = seaStateKeys{k};
    ss = str2double(key(3:end));
    fig = figure('Visible','off'); hold on
    for f = 1:numel(families)
        fam = families{f};
        if ~isfield(data.(key).families, fam), continue; end
        g = data.(key).families.(fam);
        plot(g.pressureGrid/1e6, g.aveElecPow/1e3, '-o', 'Color', colors(f,:), 'DisplayName', fam);
    end
    xlabel('High pressure (MPa)'); ylabel('Average electrical power (kW)');
    title(sprintf('Sea state %d: Hs=%.2fm, Tp=%.2fs', ss, data.(key).Hs, data.(key).Tp));
    legend('Location','best'); grid on
    set(fig,'renderer','painters');
    outFile = fullfile(plotsDir, sprintf('pressureSweep_seaState%02d.png', ss));
    saveas(fig, outFile);
    close(fig);
end
fprintf('Wrote %d pressure-sweep plot(s) to %s\n', numel(seaStateKeys), plotsDir);
end
