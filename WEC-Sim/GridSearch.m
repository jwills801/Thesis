%% Loop PI gains
KpVals = logspace(5,7,10);
% KiVals = logspace(0,4,3);
KiVals = 0;

[KpMatrix,KiMatrix] = ndgrid(KpVals,KiVals);

mcr = struct();
mcr.header = {'Kp','Ki'};
mcr.cases = [KpMatrix(:), KiMatrix(:)];
save('mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);

wecSimMCR

% Analyze results

EnergyMatrix= reshape(EnergyVector,size(KpMatrix));
figure, surf(KpMatrix,KiMatrix,-EnergyMatrix/1e6), xlabel('Kp'), ylabel('Ki'), zlabel('Energy [MJ]')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

[~,I] = min(EnergyMatrix,[],"all");
KpMatrix(I)
KiMatrix(I)
-EnergyMatrix(I)/1e6


%% Loop rectifying pressure
pressureVals = linspace(6,9,11)*1e6;

mcr = struct();
mcr.header = {'pressure'};
mcr.cases = pressureVals';
save('pressure_mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);

wecSimMCR

%% Analyze results

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
