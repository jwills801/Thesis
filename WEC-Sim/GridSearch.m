%% Loop PI gains
    % Set simu.mcrMatFile = 'gains_mcr.mat';
    % Set dt = 1e-1
    % Make sure Kp and Ki are not otherwse defined
    % Set PTO = 'Continuous PI';
    % Uncomment EHA Losses in userDefinedFunctionsMCR
% close to resonance
% KpVals = linspace(2,8,7)*1e5;
% KiVals = linspace(-3,1,5)*1e5;



% far from resonance
KpVals = linspace(2,7,6)*1e6;
KiVals = linspace(0,5,6)*1e6;

% KpVals = 4e6;
% KiVals = 2e6;


[KpMatrix,KiMatrix] = ndgrid(KpVals,KiVals);

mcr = struct();
mcr.header = {'Kp','Ki'};
mcr.cases = [KpMatrix(:), KiMatrix(:)];
save('gains_mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);
ElectricEnergyVector = NaN(size(mcr.cases,1),1);

tic
wecSimMCR
toc

% Analyze results
EnergyMatrix= reshape(EnergyVector,size(KpMatrix));
figure, surf(KpMatrix,KiMatrix,-EnergyMatrix/1e3/simu.endTime), xlabel('Kp'), ylabel('Ki'), zlabel('Average Mechanical Power [kW]')
fig = gcf; set(fig,'Color', 'white');
ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
% fileName = 'figures/continuousPI/gridSearchCloseToResonance/gridSearchCloseToResonance';
fileName = 'figures/continuousPI/gridSearchFarFromResonance/gridSearchFarFromResonance';
exportgraphics(fig,[fileName,'.pdf'],'ContentType', 'vector','Resolution', 600);
exportgraphics(fig,[fileName,'.png'],'Resolution', 600);
savefig([fileName,'.fig']);

[~,I] = min(EnergyMatrix,[],"all");
KpMatrix(I)/1e5
KiMatrix(I)/1e5
-EnergyMatrix(I)/1e3/simu.endTime

% Eha Codesign
ElectricEnergyMatrix = reshape(ElectricEnergyVector,size(KpMatrix));
figure, surf(KpMatrix,KiMatrix,-ElectricEnergyMatrix/1e3/simu.endTime), xlabel('Kp'), ylabel('Ki'), zlabel('Average Electrical Power [kW]')
fig = gcf; set(fig,'Color', 'white');
ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
% fileName = 'figures/continuousPI/EhaGridSearchCloseToResonance/EhaGridSearchCloseToResonance';
fileName = 'figures/continuousPI/EhaGridSearchFarFromResonance/EhaGridSearchFarFromResonance';
exportgraphics(fig,[fileName,'.pdf'],'ContentType', 'vector','Resolution', 600);
exportgraphics(fig,[fileName,'.png'],'Resolution', 600);
savefig([fileName,'.fig']);
[~,I] = min(ElectricEnergyMatrix,[],"all");
KpMatrix(I)/1e5
KiMatrix(I)/1e5
-ElectricEnergyMatrix(I)/1e3


%% Loop rectifying pressure
    % Set simu.mcrMatFile = 'pressure_mcr.mat';
    % Set dt = 1e-2
    % Make sure pressure is not otherwse defined
    % Set PTO = 'Rectifying';
    % Comment EHA Losses in userDefinedFunctionsMCR
pressureVals = linspace(10,20,11)*1e6;
% pressureVals = linspace(3,9,7)*1e6;

mcr = struct();
mcr.header = {'pressure'};
mcr.cases = pressureVals';
save('pressure_mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);
tic

wecSimMCR
toc

%  Analyze results
figure, plot(pressureVals/1e6,-EnergyVector/1e3/simu.endTime), xlabel('Pressure [MPa]'), ylabel('Average Mechanical Power [kW]'),grid 
fig = gcf; set(fig,'Color', 'white');
ax = findobj(fig, 'Type', 'axes'); set(ax,'FontSize', 12,'LineWidth', 2,'FontWeight', 'bold');                 % Box around axes
lines = findobj(ax, 'Type', 'line'); set(lines, 'LineWidth', 3);
% fileName = 'figures/checkValve/gridSearchCloseToResonance/gridSearchCloseToResonance';
fileName = 'figures/checkValve/gridSearchFarFromResonance/gridSearchFarFromResonance';
% exportgraphics(fig,[fileName,'.pdf'],'ContentType', 'vector','Resolution', 600);
% exportgraphics(fig,[fileName,'.png'],'Resolution', 600);
% savefig([fileName,'.fig']);

[~,I] = min(EnergyVector);
pressureVals(I)/1e6
-EnergyVector(I)/1e3
