function out = getSpectrum(BEMdata,Tp,Hs)
% Author: Jackson Wills
% Date made: Feb 23, 2026

% Description: Generate Pierson-Moskowitz Spectrum and the average wave
% resource available.

% Code dependencies:
    % CalculateWaveNumber (function defined in this file)
    % CalculateWavePower (function defined in this file)


% Define frequencies to consider
omegaMin = min(BEMdata.hydro.w); %rad/s
omegaMax = max(BEMdata.hydro.w); % rad/s
omegaNum = length(BEMdata.hydro.w); % number of frequencies considered
omegaVals = linspace(omegaMin,omegaMax,omegaNum);

% Pierson-Moskowitz Spectrum for fully develped sea states
B = 5/(4*Tp^4);
A = B*Hs^2/4;
f = omegaVals/(2*pi); % Frequency [Hz]
fSpectrum = A./(f.^5).*exp(-B./f.^4);
spectrum = fSpectrum./(2*pi); 

% output values
out.S_w = spectrum;
out.w = omegaVals;
out.dw = [0, diff(omegaVals)];


% Plot spectrum
% figure, semilogx(f,fSpectrum), grid, xlabel('Frequency [Hz]'), ylabel('Spectrum [m^2/Hz]')

end