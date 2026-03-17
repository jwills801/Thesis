function optTraj = getOptimal(params,wave)
t = params.simu.time;

% Load info about the system
sys = params.phys.sys;
F = wave.torque.elevation2torque;

% load info about the wave
S_w = wave.spectrum.S_w;
w = wave.spectrum.w;
dw = wave.spectrum.dw;
phase = wave.torque.phase;
ramp = wave.torque.ramp;

% Calculate impedance
Z = freqresp(sys,w);
Rr = real(1./squeeze(Z(1,1,:)))';

% optimal power output
avePow = trapz(w,abs(F).^2.*S_w*2/8./Rr);

% Trajectory
thetaDot = ramp.* real(sum( exp(1i*(t*w+phase)) .* ...
    (F/2./Rr.*sqrt(2*S_w.*dw)) ,2));

% For regular waves
    % H = freqresp(params.sys,2*pi/params.period);
    % thetaDotOpt = Texc/2/real(1/H(1));

% Derivative and sntegral of signal 
theta = cumtrapz(t,thetaDot);
thetaDDot = [0;diff(thetaDot)./diff(t)];

% output results
optTraj.theta = theta;
optTraj.thetaDot = thetaDot;
optTraj.thetaDDot = thetaDDot;
optTraj.time = t;
optTraj.avePow = avePow;
end