function out = getTorqueTimeSeries(params,wave)

% load parameters
g = params.phys.g;
rho = params.phys.rho;
w = wave.spectrum.w;
dw = wave.spectrum.dw;
S_w = wave.spectrum.S_w;

% Load time vector
t = params.simu.time;
rampTime = params.simu.rampTime;


% Generate ramp function
ramp = .5*(1+cos(pi + pi/rampTime*t)).*(t<rampTime) + (t>=rampTime);


% Randomly assign phases to each frequencies
rng(1) % Sets random seed for repeatability
phase = 2*pi*rand(size(w));

% Interpolate from the frequencies used in capytaine to the ones we want to
% use here
Kexc = interp1(wave.w,wave.Kexc,w);

% Clean the real data
    % If any entries are infinite, they should be zero
Kexc(~isfinite(Kexc))=0;

% Put it all together to get the excitation torque
Texc = ramp.* real(sum( exp(1i*(t*w+phase)) .* (Kexc.*sqrt(2*S_w.*dw)) ,2));

% output
out.Texc = Texc;
out.time = t;
out.elevation2torque = Kexc;
out.phase = phase;
out.ramp = ramp;

% Plot results
% figure, plot(t,Texc), ylabel('Exciting Torque [Nm]'), xlabel('Time [s]')

end