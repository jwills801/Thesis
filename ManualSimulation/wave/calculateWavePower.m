function power = calculateWavePower(params,wave)
% load constants
g = params.phys.g;
rho = params.phys.rho;
w = wave.spectrum.w;
dw = wave.spectrum.dw;
d = params.phys.waterDepth;
S_w = wave.spectrum.S_w;

% Calculate wave number
waveNumber = CalculateWaveNumber(w,d);

% Calculate Wave Power
power = sum((1/2)*rho*g*S_w.*dw.*sqrt(g./waveNumber.*tanh(waveNumber.*d)).*(1 + 2.*waveNumber.*d./sinh(2.*waveNumber.*d)));

end

function waveNumber = CalculateWaveNumber(omegaVals,waterDepth)
waveNumber = NaN(size(omegaVals));
for omegaInd = 1:length(omegaVals)
    omega = omegaVals(omegaInd);
    g = 9.81;

    % Deep water approximation, initial guess
    k = omega.^2./g;

    % Iterate for shallow and intermediate water, full dispersion relationship
    if isfinite(waterDepth)
        for i = 1:100
            k = omega.^2./g./tanh(k.*waterDepth);
        end
    end
    waveNumber(omegaInd) = k;
end
% plot
% figure, plot(omegaVals,waveNumber,omegaVals,omegaVals.^2./g)
% legend('Shallow water','Deep water'), grid
% xlabel('Frequency [Hz]'), ylabel('wave Number [1/m]')

end