%% Loop PI gains
KpVals = linspace(3,9,7)*1e6;
% KiVals = logspace(0,4,3);

KiVals = [-3e6 -1e6 -3e5 -1e5 0 1e5];

% KiVals = -1e6;
% KpVals = 4e6;

[KpMatrix,KiMatrix] = ndgrid(KpVals,KiVals);

mcr = struct();
mcr.header = {'Kp','Ki'};
mcr.cases = [KpMatrix(:), KiMatrix(:)];
save('gains_mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);
ElectricEnergyVector = NaN(size(mcr.cases,1),1);

wecSimMCR

%% Analyze results

%EnergyMatrix= reshape(EnergyVector,size(KpMatrix));
figure, surf(KpMatrix,KiMatrix,-EnergyMatrix/1e6), xlabel('Kp'), ylabel('Ki'), zlabel('Energy [MJ]')
ylim([-2e6 1e5])
fig = gcf; set(fig,'Color', 'white');
ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
exportgraphics(fig,'figures/continuousPI/gridSearchCloseToResonance/gridSearchCloseToResonance.pdf','ContentType', 'vector','Resolution', 600);
exportgraphics(fig,'figures/continuousPI/gridSearchCloseToResonance/gridSearchCloseToResonance.png','Resolution', 600);
savefig('figures/checkValve/gridSearchCloseToResonance/gridSearchCloseToResonance.fig');


[~,I] = min(EnergyMatrix,[],"all");
KpMatrix(I)/1e6
KiMatrix(I)/1e6
-EnergyMatrix(I)/1e6

%% analyze results for one dof
figure, plot(KiVals/1e6,-EnergyVector/1e6/40), xlabel('Ki'), ylabel('Average Mechanical Power [MW]'),grid 
% fig = gcf; set(fig,'Color', 'white');
% ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
% lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
% exportgraphics(fig,'figures/continuous/gridSearchCloseToResonance.pdf','ContentType', 'vector','Resolution', 600);
% exportgraphics(fig,'figures/checkValve/gridSearchCloseToResonance.png','Resolution', 600);
% savefig('figures/checkValve/gridSearchCloseToResonance.fig');

[~,I] = min(EnergyVector,[],"all");
KiMatrix(I)
-EnergyVector(I)/1e6


%% Loop rectifying pressure
pressureVals = linspace(6,9,11)*1e6;

mcr = struct();
mcr.header = {'pressure'};
mcr.cases = pressureVals';
save('pressure_mcr.mat','mcr')

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
