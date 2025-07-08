% Run all PTO cases
% PTO_wave_mcr.mat
% simu.dt = 1e-3;
% Uncomment peakWavePeriod, codesign, and PTO

clear, close all
PTO_options = {'Active Valving', 'HHEA', 'Passive Valving','EHA'};
peakWavePeriod_options = [8, 20]; 
codesign_options = 0;

% Create combinations (4 PTOs × 2 periods = 8 cases)
[PTO_grid, period_grid] = meshgrid(PTO_options, peakWavePeriod_options);
combinations = [PTO_grid(:), num2cell(period_grid(:)),num2cell(zeros(size(PTO_grid(:))))];  % Convert to cell array

% Add EHA cases with codesign = 1
eha_idx = strcmp(PTO_grid(:), 'EHA');
eha_cases = [PTO_grid(eha_idx), num2cell(period_grid(eha_idx)), num2cell(ones(sum(eha_idx), 1))];
combinations = [combinations; eha_cases];

mcr.header = {'PTO', 'peakWavePeriod','codesign'};  % Column names
mcr.cases = combinations;  % Each row is a test case
save('PTO_wave_mcr.mat', 'mcr');

EnergyVector = NaN(size(mcr.cases,1),1);
outEnergyVector = NaN(size(mcr.cases,1),1);

%% run all cases
for imcr = 1:length(mcr.cases)
    for variableNameInd = 1:length(mcr.header)
        eval([mcr.header{variableNameInd} ,'= mcr.cases{imcr,variableNameInd}'])
    end
    wecSim
    EnergyVector(imcr) = Energy;
    outEnergyVector(imcr) = outputEnergy;
end
