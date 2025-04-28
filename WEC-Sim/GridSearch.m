%% Loop PI gains
KpVals = logspace(6,8,10);
KiVals = logspace(6,8,10);

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
pressureVals = linspace(4,10,5)*1e6;

mcr = struct();
mcr.header = {'pressure'};
mcr.cases = pressureVals';
save('pressure_mcr.mat','mcr')

EnergyVector = NaN(size(mcr.cases,1),1);

wecSimMCR

% Analyze results

figure, plot(pressureVals/1e6,-EnergyVector/1e6), xlabel('Pressure [MPa]'), ylabel('MJ')

[~,I] = min(EnergyVector);
pressureVals(I)/1e6
-EnergyVector(I)/1e6
