hydro = struct();
hydro = readCAPYTAINE(hydro,'System_Identification/Capytaine/oswec_new2.nc');

params = struct();
params.hydro = hydro;
params.rho = 1025;
params.g = 9.81;

params.Ts = 20;
params.Hs = 0.99;
params.alpha = NaN;
params.rampTime = 20;
params.omega = waves.omega'; params.dw = waves.dOmega';
params.phase = waves.phase';
[params.spectrum,power] = JS_Spectrum(params.Ts,params.Hs,params.alpha,params.omega);
    % figure, plot(omega,spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)

%% Simulate Dynamics (euler integration)
params.dt = .01;
params.finalTime = 10;
params.time = 0:dt:finalTime;

position = NaN(size(time));
velocity = NaN(size(time));
position(1) = 0;
velocity(1) = 0;

timeInd = 1;

out = WECDynamics(position, velocity,params,time(timeInd),timeInd)



%% Functions

function out = WECDynamics(position, velocity, obj,time,timeInd)

% Excitation Torque
torqueExcitation = GetExcitationForce(obj,time);

% Stiffness Torque
torqueStiffness = position*obj.hydro.Khs*obj.rho*obj.g;

% Save to output
out.torqueExcitation = torqueExcitation;
end

function torqueExcitation = GetExcitationForce(obj,time)
w_vals = obj.hydro.w;
ex_re_vals = squeeze(obj.hydro.ex_re(5,1,:));
ex_im_vals = squeeze(obj.hydro.ex_im(5,1,:));
omega = obj.omega;
ex_re = interp1(w_vals,ex_re_vals,omega);
ex_im = interp1(w_vals,ex_im_vals,omega);

F = squeeze((ex_re+ex_im*1i)) * obj.rho*obj.g;
ramp = .5*(1+cos(pi + pi/obj.rampTime*time)).*(time<obj.rampTime) + (time>=obj.rampTime);
torqueExcitation = ramp.* real(sum( exp(1i*(time*omega+obj.phase)) .* (F.*sqrt(2*obj.spectrum.*obj.dw)) ,2));
end

function [spectrum,power] = JS_Spectrum(Tp,Hs,gamma,omega)
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