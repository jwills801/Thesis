hydro = struct();
hydro = readCAPYTAINE(hydro,'System_Identification/Capytaine/oswec_new2.nc');

params = struct();
params.hydro = hydro;
params.rho = 1025;
params.g = 9.81;

centerOfGravity = [0;0;-3.9;0;0;0];
centerOfBuoyancy = [0;0;-4.5915;0;0;0];

position = output.bodies(1).position';
velocity = output.bodies(1).velocity';
acceleration = output.bodies(1).acceleration';

forceTotal = output.bodies(1).forceTotal';
forceExcitation = output.bodies(1).forceExcitation';
forceAddedMass = output.bodies(1).forceAddedMass';
forceLinearDamping = output.bodies(1).forceLinearDamping';
forceMorisonAndViscous = output.bodies(1).forceMorisonAndViscous';
forceRadiationDamping= output.bodies(1).forceRadiationDamping';
forceRestoring = output.bodies(1).forceRestoring';
forceTotalCheck = forceExcitation - forceAddedMass - forceLinearDamping - forceMorisonAndViscous - forceRadiationDamping - forceRestoring;
% figure, plot(time,forceTotal-forceTotalCheck)

params.Ts = 20;
params.Hs = 0.99;
params.alpha = NaN;
params.rampTime = 20;
params.omega = waves.omega'; params.dw = waves.dOmega';
params.phase = waves.phase';
[params.spectrum,power] = JS_Spectrum(params.Ts,params.Hs,params.alpha,params.omega);
    % figure, plot(omega,spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)

%%
%% Radiation Damping
rad = struct();
rad.time = 0:.01:50;
rad.w = linspace(min(hydro.w)*.01,max(hydro.w),1001);

B = interp1(hydro.w,squeeze(hydro.B(5,5,:)).*hydro.w'*params.rho,rad.w); % This extra w and rho is to dimensionalize B 
B(~isfinite(B)) = 0;
rad.Kr = (2/pi)*trapz(rad.w,B.*(cos(rad.w.*rad.time(:))), 2);    

%%
% Convolution Integral
rad.convTime = 0:.01:30; % This time matches the time wecsim uses
rad.convolutionLength = 30; % It takes about 30 seconds for Kr terms to decay
F = NaN(length(time),1);   % Initilize force
timeInds = 3401:3501; timeInds = 2:1001;
figWaitbar = waitbar(0,'Convolution Integral');
for time_ind = timeInds%:length(time)         % Loop over time
    waitbar((time_ind-timeInds(1))/(timeInds(end)-timeInds(1)),figWaitbar)
        % Find the start and ending indices of the fine radiation time vector
    if time(time_ind) > rad.convolutionLength
        startTime = time(time_ind)-rad.convolutionLength; % Start time of convolution integral
        [~,  endInd] = min(abs(rad.convTime-rad.convolutionLength)); % ending index corresponding to rad.convolutionLength seconds of integration
    else
        startTime = 0;
        [~,  endInd] = min(abs(rad.convTime-time(time_ind))); % ending index corresponding to length of integration
    end

        % Compute integrand for each value of tau
    integrand = NaN(endInd,1);         % Initialize integrand
    for tau_ind = 1:endInd              % Loop over values of tau
        tau = startTime+rad.convTime(tau_ind);
        v = interp1(time,velocity(5,:),tau); % interpolate velocity since it is on a coarser time mesh
        Kr = interp1(rad.time,rad.Kr,time(time_ind)-tau);
        integrand(tau_ind) = Kr*v;           % Store the integrand for each value of tau
    end

    % integrate over tau
    F(time_ind) = trapz(rad.convTime(1:endInd),integrand); % perform the integration
end
close(figWaitbar)
r = 5; theta = position(5,:); forces = forceRadiationDamping;
torqueRadiationDamping = forces(5,:)  + r*forces(1,:).*cos(theta) - r*forces(3,:).*sin(theta);
figure, plot(time,F,time,torqueRadiationDamping), xlim([time(timeInds(1))-5, time(timeInds(end))+5])

%% Radiation Damping State Space
figure, plot(rad.time,rad.Kr), hold on, ax1 = gca; fig1 = gcf; grid
figure, plot(time,F), hold on, ax2 = gca; grid
for n = 3
    theta_0 = ones(2*n-1,1);
    fun = @(theta) CompareImpulse(theta,rad,n);
    options = optimoptions('lsqnonlin','PlotFcn','optimplotresnorm');
    theta = lsqnonlin(fun,theta_0,[],[],options)

    [error,K_hat] = CompareImpulse(theta,rad,n);
    plot(ax1,rad.time,K_hat)

    A = [zeros(1,n-1), -theta(1);
        eye(n-1), -theta(2:n)];
    B = [0;theta(n+1:end)]*1e7;
    C = [zeros(1,n-1), 1];
    rad.ss = ss(A,B,C,0);
    forceRadiationDampingCheck = lsim(rad.ss,velocity(5,:),time);
    plot(ax2,time,forceRadiationDampingCheck)
end
hold(ax1,'off'), xlabel(ax1,'Time [s]'), ylabel(ax1,'Radiation Impulse Response [Nm/rad]'), xlim([0 20])
% exportgraphics(fig1, 'figures/StateSpaceRadDamping.pdf', 'ContentType', 'vector');
% savefig(fig1,"figures/StateSpaceRadDamping.fig")

hold(ax2,'off')

[mag, phase] = bode(rad.ss,w);
figure, subplot(211), semilogx(w,20*log10(abs(K)),w,20*log10(abs(squeeze(mag)))), grid, xlabel('Frequency [rad/s]'), ylabel('Magnitude [dB]'), xlim([.1 10])
        subplot(212), semilogx(w,180/pi*angle(K),w,squeeze(phase)), grid, xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'), xlim([.1 10])
% exportgraphics(gcf, 'figures/StateSpaceRadDamp_bodeCheck.pdf', 'ContentType', 'vector');
% savefig("figures/StateSpaceRadDamp_bodeCheck.fig")

figure, plot(time, torqueRadiationDamping,time,F,time,forceRadiationDampingCheck)
figure, bode(rad.ss), hold on, bode(Khat), hold off
return
%%
w = hydro.w';
B = squeeze(hydro.B(5,5,:)).*w*params.rho; B(~isfinite(B)) = 0;
A = squeeze(hydro.A(5,5,:))*params.rho;
Ainf = hydro.Ainf(5,5)*params.rho;

Kr = (2/pi)*trapz(w,B.*(cos(w.*rad.time)),1)';
A2 = NaN(size(w));
for wInd = 1:length(w)
    A2(wInd) = Ainf - 1/w(wInd)* trapz(rad.time,Kr.*sin(w(wInd)*rad.time'));
end
% figure, plot(hydro.w,A,hydro.w,A2)

K = B + 1j*w.*(A-Ainf);
K = B + 1j*w.*(A2-Ainf);

w2 = logspace(-1,1,100);
K2 = interp1(w,K,w2);
figure, plot(w,K,w2,K2)

% figure, semilogx(w,20*log10(abs(K)),w,20*log10(abs(K2)))

n = 1; m = 2;

[b,a] = invfreqs(K,w,n,m);
%b(2) = 0;
sys = tf(b,a);

[b2,a2] = invfreqs(K2,w2,n,m);
b2(2) = 0;
sys2 = tf(b2,a2);

[mag,phase] = bode(sys,w);
K_est = squeeze(mag);
% figure, semilogx(w,20*log10(abs(K)),w,20*log10(abs(K_est)))
[mag2,phase2] = bode(sys2,w);
K_est2 = squeeze(mag2);
figure, semilogx(w2,20*log10(abs(K2)),w,20*log10(abs(K_est2)))

[y,tOut]= impulse(sys);
[y2,tOut2]= impulse(sys2);
% figure, plot(rad.time,rad.Kr,tOut,y,tOut2,y2)
%%
%
theta_0 = ones(3,1);
func2 = @(theta) fun2(theta,K2,w2);
options = optimoptions('lsqnonlin','PlotFcn','optimplotresnorm','FunctionTolerance',1e-6);
theta = lsqnonlin(func2,theta_0,[],[],options);
% theta = [6;2.3; 6.8];
s = 1j*w;
K_hat = theta(1)*1e8*s ./ (s.^2 + theta(2)*s + theta(3) );
%
% theta = [0.5;0.3;1.25];
K_hat = theta(1)*1e8*s ./ (s.^2 + 2*theta(2)*theta(3)*s + theta(3)^2 );
figure, subplot(211), semilogx(w,20*log10(abs(K)),w,20*log10(abs(K_hat))), grid, xlabel('Frequency [rad/s]'), ylabel('Magnitude [dB]'), xlim([.1 10])
        subplot(212), semilogx(w,180/pi*angle(K),w,180/pi*angle(K_hat)), grid, xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'), xlim([.1 10])
% exportgraphics(gcf, 'figures/RadDampTF.pdf', 'ContentType', 'vector');
% savefig("figures/RadDampTF.fig")

sys = tf([theta(1)*1e8,0],[1,theta(2:3)']);
sys = tf([theta(1)*1e8,0],[1,2*theta(2)*theta(3),theta(3)^2])
[y,tOut]= impulse(sys,rad.time);
figure, plot(rad.time,rad.Kr,tOut,y), grid, xlabel('Time [s]'), ylabel('Radiation Impulse Response [Nm/rad]'), xlim([0 20])
% exportgraphics(gcf, 'figures/RadDampingTf_impulseCheck.pdf', 'ContentType', 'vector');
% savefig("figures/StateSpaceRadDampingTf_impulseCheck.fig")

forceRadiationDampingCheckTf = lsim(sys,velocity(5,:),time);
figure, plot(time,F,time,forceRadiationDampingCheck,time,forceRadiationDampingCheckTf), grid
legend('Convolution','Time Domain','Frequency Domain'),xlim([time(timeInds(1)), time(timeInds(end))])
ylabel('Radiation Torque [Nm]'), xlabel('Time [s]')

%%
ex_re = squeeze(hydro.ex_re(5,1,:));
ex_im = squeeze(hydro.ex_im(5,1,:));
F = (ex_re+ex_im*1i) * rho*g;
I = 5.025e6;
Khs = (rho*Vol*4.5-m*5)*g;
elevation2magnitude = F ./ (-w.^2.*(I+A2) + 1i*w.*B + Khs);
figure, subplot(211), semilogx(w,20*log10(abs(elevation2magnitude))), grid, xlabel('Frequency [rad/s]'), ylabel('Magnitude [dB]'), xlim([.1 10])
        subplot(212), semilogx(w,180/pi*angle(elevation2magnitude)), grid, xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'), xlim([.1 10])

torque2velocity = 1 ./ (-w.^2.*(I+A2) + 1i*w.*B + Khs);
figure, subplot(211), semilogx(w,20*log10(abs(torque2velocity))), grid, xlabel('Frequency [rad/s]'), ylabel('Magnitude [dB]'), xlim([.1 10])
        subplot(212), semilogx(w,180/pi*angle(torque2velocity)), grid, xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'), xlim([.1 10])

return

%% RM3 Frequency Response
%RM3 = struct();
%RM3 = readWAMIT(RM3,'WEC-Sim-source/examples/RM3/hydroData/rm3.out',[]);
RM3.Bd = squeeze(RM3.B(3,3,:)).*RM3.w'*rho;
RM3.Ad = squeeze(RM3.A(3,3,:))*rho;
RM3.K = RM3.Khs(3,3,1)*rho*g;
RM3.m = 727e3;
RM3.torque2velocity = 1 ./ (-RM3.w'.^2.*(RM3.m+RM3.Ad) + 1i*RM3.w'.*RM3.Bd + RM3.K);
figure, semilogx(w,20*log10(abs(torque2velocity)),RM3.w,20*log10(abs(RM3.torque2velocity))), grid, xlabel('Frequency [rad/s]'), ylabel('Magnitude [dB]'), xlim([.1 5]), legend('OSWEC','RM3')
        %subplot(212), semilogx(w,180/pi*angle(torque2velocity),RM3.w,180/pi*angle(RM3.torque2velocity)), grid, xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'), xlim([.1 10])

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

function [error,K_hat] = CompareImpulse(theta,rad,n)
    A = [zeros(1,n-1), -theta(1);
        eye(n-1), -theta(2:n)];
    B = [0;theta(n+1:end)]*1e7;
    C = [zeros(1,n-1), 1];

K_hat = NaN(length(rad.time),1);
for time_ind = 1:length(rad.time)
    K_hat(time_ind) = C*expm(A*rad.time(time_ind))*B;
end
error = rad.Kr - K_hat;
end

function error = fun2(theta,K,w)
s = 1j*w;
K_hat = theta(1)*1e8*s ./ (s.^2 + theta(2)*s + theta(3) );
K_hat = theta(1)*1e8*s ./ (s.^2 + 2*theta(2)*theta(3)*s + theta(3)^2 );
error = abs(K - K_hat);
end