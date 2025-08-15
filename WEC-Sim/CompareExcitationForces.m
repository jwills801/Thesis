hydro = struct();
hydro = readCAPYTAINE(hydro,'System_Identification/Capytaine/oswec_new2.nc');

rho = 1025;
g = 9.81;

Ts = 20;
Hs = 0.99;
alpha = NaN;
rampTime = 20;
omega = waves.omega'; dw = waves.dOmega';
phase = waves.phase';
[spectrum,power] = JS_Spectrum(Ts,Hs,alpha,omega);
    % figure, plot(omega,spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)

ex_re = interp1(hydro.w,squeeze(hydro.ex_re(5,1,:)),omega);
ex_im = interp1(hydro.w,squeeze(hydro.ex_im(5,1,:)),omega);
F = squeeze((ex_re+ex_im*1i)) * rho*g;
ramp = .5*(1+cos(pi + pi/rampTime*time)).*(time<rampTime) + (time>=rampTime);
torqueExcitation = ramp.* real(sum( exp(1i*(time*omega+phase)) .* (F.*sqrt(2*spectrum.*dw)) ,2));
figure, plot(time,torqueExcitation)

% 6 DOF version of excitaton torque
hydro6DOF = struct();
hydro6DOF = readCAPYTAINE(hydro6DOF,'System_Identification/Capytaine/oswec_new.nc');
hydro6DOF.forceExcitation = NaN(6,length(time));
for i = 1:sum(hydro6DOF.dof)
    ex_re = interp1(hydro6DOF.w,squeeze(hydro6DOF.ex_re(i,1,:)),omega);
    ex_im = interp1(hydro6DOF.w,squeeze(hydro6DOF.ex_im(i,1,:)),omega);
    hydro6DOF.F = squeeze((ex_re'+ex_im'*1i)) * rho*g;
    hydro6DOF.forceExcitation(i,:) = ramp.*real(exp(1i*(time*omega+phase)) * (hydro6DOF.F.*sqrt(2*spectrum'.*dw')));
end
figure, % Compares these to wecsim (which uses oswec_full.nc)
    subplot(131), plot(time,forceExcitation(1,:),time,hydro6DOF.forceExcitation(1,:)), title('Surge'), grid, ylabel('Force or Torque [N or Nm]')
    subplot(132), plot(time,forceExcitation(3,:),time,hydro6DOF.forceExcitation(3,:)), title('Heave'), grid
    subplot(133), plot(time,forceExcitation(5,:),time,hydro6DOF.forceExcitation(5,:)), title('Pitch'), grid, xlabel('Time [s]')
theta = 0;
hydro6DOF.torqueExcitation = hydro6DOF.forceExcitation(5,:)  + r*hydro6DOF.forceExcitation(1,:).*cos(theta) - r*hydro6DOF.forceExcitation(3,:).*sin(theta);
figure, plot(time,hydro6DOF.torqueExcitation,time,torqueExcitation) % It matches!!







function [spectrum,power] = JS_Spectrum(Tp,Hs,gamma,omega)
% Inputs:
    % omega: 1 x nOmega frequencies in rad/s
waterDepth = 10.9;
g = 9.81; rho = 1023;
frequency = omega'/(2*pi);
dw = [omega(2)-omega(1), diff(omega)];

% Pierson-Moskowitz Spectrum from IEC TS 62600-2 ED2 Annex C.2 (2019)
bPM = (5/4)*(1/Tp)^(4);
aPM =  bPM*(Hs/2)^2;
fSpectrum  = (aPM*frequency.^(-5).*exp(-bPM*frequency.^(-4)));            % Wave Spectrum [m^2-s] for 'EqualEnergy'

% JONSWAP Spectrum from IEC TS 62600-2 ED2 Annex C.2 (2019)
fp = 1/Tp;
siga = 0.07;sigb = 0.09;                                    % cutoff frequencies for gamma function
[lind,~] = find(frequency<=fp);
[hind,~] = find(frequency>fp);
gammaAlpha = zeros(size(frequency));
gammaAlpha(lind) = gamma.^exp(-(frequency(lind)-fp).^2/(2*siga^2*fp^2));
gammaAlpha(hind) = gamma.^exp(-(frequency(hind)-fp).^2/(2*sigb^2*fp^2));
C = 1 - 0.287*log(gamma);
%fSpectrum = C*fSpectrum.*gammaAlpha;                                % Wave Spectrum [m^2-s]
spectrum = fSpectrum'/(2*pi);

waveNumber = calcWaveNumber(omega,waterDepth,g,0); %Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
power = sum((1/2)*rho*g*spectrum.*dw.*sqrt(g./waveNumber.*tanh(waveNumber.*waterDepth)).*(1 + 2.*waveNumber.*waterDepth./sinh(2.*waveNumber.*waterDepth)));

end