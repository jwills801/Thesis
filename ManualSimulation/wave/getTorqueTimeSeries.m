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
ex_re = interp1(wave.BEMdata.hydro.w,squeeze(wave.BEMdata.hydro.ex_re(5,1,:)),w);
ex_im = interp1(wave.BEMdata.hydro.w,squeeze(wave.BEMdata.hydro.ex_im(5,1,:)),w);

% Clean the real data
    % the first few entries in this are infinite, but they should be zero
ex_re(~isfinite(ex_re))=0;

% Put together the real and imaginary components and dimensionalize
F = (ex_re+ex_im*1i) * rho*g;

% Put it all together to get the excitation torque
Texc = ramp.* real(sum( exp(1i*(t*w+phase)) .* (F.*sqrt(2*S_w.*dw)) ,2));

% output
out.Texc = Texc;
out.time = t;
out.elevation2torque = F;
out.phase = phase;
out.ramp = ramp;

% Plot results
% figure, plot(t,Texc), ylabel('Exciting Torque [Nm]'), xlabel('Time [s]')

end