clear, close all
PTO_options = {'Active Valving', 'EHA', 'HHEA', 'Passive Valving'};
peakWavePeriod_options = [8, 20]; 
peakWavePeriod_options = [20]; 

% Create combinations (4 PTOs × 2 periods = 8 cases)
[PTO_grid, period_grid] = meshgrid(PTO_options, peakWavePeriod_options);
combinations = [PTO_grid(:), num2cell(period_grid(:))];  % Convert to cell array

% Build the mcr struct
mcr = struct();
mcr.header = {'PTO', 'peakWavePeriod'};  % Column names
mcr.cases = combinations;  % Each row is a test case

% Save to .mat file
save('PTO_wave_mcr.mat', 'mcr');

EnergyVector = NaN(size(mcr.cases,1),1);

wecSimMCR

% Analyze results

figure, plot(pressureVals/1e6,-EnergyVector/1e6/40), xlabel('Pressure [MPa]'), ylabel('Average Mechanical Power [MW]'),grid 
fig = gcf; set(fig,'Color', 'white');
ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
exportgraphics(fig,'figures/checkValve/gridSearchCloseToResonance.pdf','ContentType', 'vector','Resolution', 600);
exportgraphics(fig,'figures/checkValve/gridSearchCloseToResonance.png','Resolution', 600);
savefig('figures/checkValve/gridSearchCloseToResonance.fig');

[~,I] = min(EnergyVector);
pressureVals(I)/1e6
-EnergyVector(I)/1e6