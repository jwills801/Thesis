% The goal of this code to recreate the dynamics in WECSim

% Data from WEC-Sim
r = 5; % Distance from hinge to center of gravity
time = output.bodies(1).time;
rho = 1000;
g = 9.81;
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

%% Data from Capytaine
hydro = struct();
hydro = readCAPYTAINE(hydro,'System_Identification/Capytaine/oswec_full.nc');
hydro = radiationIRF(hydro,[],[],[],[],[]);

hydro_new = struct();
hydro_new = readCAPYTAINE(hydro_new,'System_Identification/Capytaine/oswec_new.nc');

hydro_new2 = struct();
hydro_new2 = readCAPYTAINE(hydro_new2,'System_Identification/Capytaine/oswec_new2.nc');
figure,
    subplot(131), plot(hydro.w,squeeze(hydro.B(1,1,:)),hydro_new.w,squeeze(hydro_new.B(1,1,:)),hydro_new2.w,squeeze(hydro_new2.B(1,1,:))), title('Surge'), grid, ylabel('Damping Coeff')
    subplot(132), plot(hydro.w,squeeze(hydro.B(3,3,:)),hydro_new.w,squeeze(hydro_new.B(3,3,:)),hydro_new2.w,squeeze(hydro_new2.B(3,3,:))), title('Heave'), grid
    subplot(133), plot(hydro.w,squeeze(hydro.B(5,5,:)),hydro_new.w,squeeze(hydro_new.B(5,5,:)),hydro_new2.w,squeeze(hydro_new2.B(5,5,:))), title('Pitch'), grid, xlabel('Frequency [rad/s]')

figure,
    subplot(231), plot(hydro.w,squeeze(hydro.ex_re(1,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(1,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_re(1,1,:))), title('Surge'), grid, ylabel('Real Exitation Component')
    subplot(232), plot(hydro.w,squeeze(hydro.ex_re(3,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(3,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_re(3,1,:))), title('Heave'), grid
    subplot(233), plot(hydro.w,squeeze(hydro.ex_re(4,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(5,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_re(5,1,:))), title('Pitch'), grid, xlabel('Frequency [rad/s]')
    subplot(234), plot(hydro.w,squeeze(hydro.ex_im(1,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(1,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_im(1,1,:))), title('Surge'), grid, ylabel('Imaginary Exitation Component')
    subplot(235), plot(hydro.w,squeeze(hydro.ex_im(3,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(3,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_im(3,1,:))), title('Heave'), grid
    subplot(236), plot(hydro.w,squeeze(hydro.ex_im(4,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(5,1,:)),hydro_new2.w,squeeze(hydro_new2.ex_im(5,1,:))), title('Pitch'), grid, xlabel('Frequency [rad/s]')

figure,
    subplot(231), plot(hydro.w,squeeze(hydro.ex_re(1,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(1,1,:))), title('Surge'), grid, xlabel('Frequency [rad/s]') , ylabel('Real Exitation Component')
    subplot(232), plot(hydro.w,squeeze(hydro.ex_re(3,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(3,1,:))), title('Heave'), grid, xlabel('Frequency [rad/s]')
    subplot(233), plot(hydro.w,squeeze(hydro.ex_re(4,1,:)),hydro_new.w,squeeze(hydro_new.ex_re(5,1,:))), title('Pitch'), grid, xlabel('Frequency [rad/s]')
    subplot(234), plot(hydro.w,squeeze(hydro.ex_im(1,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(1,1,:))), title('Surge'), grid, xlabel('Frequency [rad/s]'), ylabel('Imaginary Exitation Component')
    subplot(235), plot(hydro.w,squeeze(hydro.ex_im(3,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(3,1,:))), title('Heave'), grid, xlabel('Frequency [rad/s]')
    subplot(236), plot(hydro.w,squeeze(hydro.ex_im(4,1,:)),hydro_new.w,squeeze(hydro_new.ex_im(5,1,:))), title('Pitch'), grid, xlabel('Frequency [rad/s]')

%% Added Mass
hydro.forceAddedMass = hydro.A(1:6,1:6,end)*acceleration*rho;
plotForces(time,forceAddedMass,hydro.forceAddedMass)

%% Hydrostatic Stiffness
hydro.forceRestoring = hydro.Khs(:,:,1)*rho*g*(position-centerOfGravity) + [0;0;g*127000-rho*g*body(1).volume;0;0;0];
plotForces(time,forceRestoring,hydro.forceRestoring)

%% Excitation Force
% phaseSeed = 1; rng(phaseSeed);
% phase = 2*pi*rand(1,length(hydro.w));
phase = waves.phase';
omega = waves.omega';
[spectrum,power] = JS_Spectrum(20,.99,3.3,omega);
dw = [omega(2)-omega(1); diff(omega')];
% figure, plot(omega,spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)
% plotSpectrum(waves)

hydro.forceExcitation = NaN(size(hydro.forceRestoring));
for i = 1:sum(hydro.dof)
    ex_re = interp1(hydro.w,squeeze(hydro.ex_re(i,1,:)),omega);
    ex_im = interp1(hydro.w,squeeze(hydro.ex_im(i,1,:)),omega);
    F = squeeze((ex_re'+ex_im'*1i)) * rho*g;
    hydro.forceExcitation(i,:) = real(exp(1i*(time*omega+phase)) * (F.*sqrt(2*spectrum'.*dw)));
end
plotForces(time,forceExcitation,hydro.forceExcitation)

[waves.power,power]

%% Radiation Damping
rad = struct();
rad.time = 0:.01:50;
rad.w = linspace(min(hydro.w),max(hydro.w),1001);
rad.Kr = NaN(sum(hydro.dof), sum(hydro.dof), length(rad.time));
rad.Kr2 = NaN(sum(hydro.dof), sum(hydro.dof), length(rad.time));
for i = 1:sum(hydro.dof)
    for j = 1:sum(hydro.dof)
        B = interp1(hydro.w,squeeze(hydro.B(i,j,:)),rad.w);
        rad.Kr(i,j,:) = (2/pi)*trapz(rad.w,B.*(cos(rad.w.*rad.time(:)).*rad.w), 2);    % This extra w is to dimensionalize B (we still need to multiply by density
    end
end
figure,
dof = 1; subplot(131), plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:))), ylabel('Surge Diagonal Entry in Kr [N/m]'), xlabel('Time')
dof = 3; subplot(132), plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:))), ylabel('Heave Diagonal Entry in Kr [N/m]'), xlabel('Time')
dof = 5; subplot(133), plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:))), ylabel('Pitch Diagonal Entry in Kr [Nm/rad]'), xlabel('Time')

% Convolution Integral
rad.convTime = 0:.1:30; % This time matches the time wecsim uses
rad.convolutionLength = 10; % It takes about 30 seconds for Kr terms to decay
F = NaN(12,length(time));   % Initilize force
for time_ind = 201:231%length(time)         % Loop over time
        % Find the start and ending indices of the fine radiation time vector
    if time(time_ind) > rad.convolutionLength
        startTime = time(time_ind)-rad.convolutionLength; % Start time of convolution integral
        [~,  endInd] = min(abs(rad.convTime-rad.convolutionLength)); % ending index corresponding to rad.convolutionLength seconds of integration
    else
        startTime = 0;
        [~,  endInd] = min(abs(rad.convTime-time(time_ind))); % ending index corresponding to length of integration
    end
        
        % Compute integrand for each value of tau
    integrand = NaN(12,endInd);         % Initialize integrand
    for tau_ind = 1:endInd              % Loop over values of tau
        tau = startTime+rad.convTime(tau_ind);
        v = interp1(time,[velocity',zeros(size(velocity'))],tau)'; % interpolate velocity since it is on a coarser time mesh
        Kr = interpolateMatrix(rad,time(time_ind)-tau)*rho;    % times by rho to dimensionalize
        % Kr = rad.Kr(:,:,endInd-tau_ind+1) * rho;          % times by rho to dimensionalize
        integrand(:,tau_ind) = Kr*v;           % Store the integrand for each value of tau
    end
        
        % integrate over tau
    F(:,time_ind) = trapz(rad.convTime(1:endInd),integrand,2); % perform the integration
end
figure, plot(time,F(5,:),time,forceRadiationDamping(5,:))

%% Inertial forces
theta = position(5,:);
torqueTotal = getNetTorque(forceTotal+forceAddedMass,theta);
mass = [diag(repmat(body(1).mass(1),1,3)), zeros(3,3);zeros(3,3), diag(repmat(body(1).inertia(1),1,3))];
I = mass+body(1).hydroData.hydro_coeffs.added_mass.inf_freq(1:6,1:6)*rho;
forceStructuralInertia = I*acceleration;
torqueStructuralInertia = getNetTorque(forceStructuralInertia,theta);
torqueStructuralInertiaLinear = (I(5,5)+I(1,1)*5^2+(I(1,5)+I(5,1))*5)*acceleration(5,:);
figure, plot(time,torqueTotal,time,torqueStructuralInertia,time,torqueStructuralInertiaLinear)

%% Recreate dynamics from scratch
copy = struct();            % structure in which all dynamics from wecsim will be copied by re calculating dynamics from scratch
copy.Hs = 0.99;             % [m] Significant Wave Height
copy.Ts = 20;               % [s] Significant Wave Period
copy.alpha = 1.7;           % Peak enhancement factor for JS spectrum
copy.time = 0:.01:time(end);           % [s] Time vector
copy.phase = waves.phase';   % use phase from wecSim
copy.omega = waves.omega';   % use frequencies from wecSim
copy.rampTime = 20;
copy.dof = 6;
copy.rho = 1023;
copy.g = 9.81;
copy.cog = [0;0;-3.9;0;0;0];
copy.mass = 127000;
copy.inertia = 1.85e6;
copy.volume = 297.3760;
copy.convolutionLength = 10;
copy.convolutionDt = .01; % This is time that I used in wecsim, but it needs to be decreased

% Initialize arrays
copy.position = NaN(copy.dof,length(copy.time));
copy.velocity = NaN(copy.dof,length(copy.time));
copy.acceleration = NaN(copy.dof,length(copy.time));

% inital conditions
copy.position(:,1) = 0;
copy.velocity(:,1) = 0;

% Calculate excitation force ahead of time
[copy.spectrum,copy.power] = JS_Spectrum(copy.Ts,copy.Hs,copy.alpha,copy.omega);
dw = [copy.omega(2)-copy.omega(1); diff(copy.omega')];
copy.dw = dw';
    % figure, plot(copy.w,copy.spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)
copy.forceExcitation = NaN(copy.dof,length(copy.time));
for i = 1:sum(copy.dof)
    ex_re = interp1(hydro.w,squeeze(hydro.ex_re(i,1,:)),copy.omega);
    ex_im = interp1(hydro.w,squeeze(hydro.ex_im(i,1,:)),copy.omega);
    copy.F(i,:) = squeeze((ex_re+ex_im*1i)) * rho*g;
end

% Calculate radiation impulse response function ahead of time 
rad = struct();
rad.time = 0:.01:copy.convolutionLength;
rad.w = linspace(min(hydro.w),max(hydro.w),1001);
rad.Kr = NaN(sum(hydro.dof), sum(hydro.dof), length(rad.time));
rad.Kr2 = NaN(sum(hydro.dof), sum(hydro.dof), length(rad.time));
for i = 1:copy.dof
    disp(['Radiation impulse response function is ',num2str(i/copy.dof*100),' % done'])
    for j = 1:copy.dof
        B = interp1(hydro.w,squeeze(hydro.B(i,j,:)),rad.w);
        rad.Kr(i,j,:) = (2/pi)*trapz(rad.w,B.*(cos(rad.w.*rad.time(:)).*rad.w), 2);    % This extra w is to dimensionalize B (we still need to multiply by density
    end
end
copy.rad = rad;

%%
copy.forceExcitation2=forceExcitation(1:6,:);
copy.forceRadiation2=forceRadiationDamping(1:6,:);
copy.forceRestoring2=forceRestoring(1:6,:);
endInd = 101;
for timeInd = 1:endInd% length(copy.time)
    timeInd
    states = [copy.position(:,timeInd);copy.velocity(:,timeInd)];
    out = WEC_dynamics(states,copy.time(timeInd),copy,hydro);
    copy.position(:,timeInd+1)  = copy.position(:,timeInd) + copy.velocity(:,timeInd)*copy.time(2);
    copy.velocity(:,timeInd+1)  = copy.velocity(:,timeInd) + out.acceleration*copy.time(2);

    % Keep outputs
    copy.forceExcitation(:,timeInd) = out.forceExcitation;
    copy.forceRadiation(:,timeInd) = out.forceRadiation;
    copy.forceRestoring(:,timeInd) = out.forceRestoring;
end
% figure, plot(copy.time,copy.position)
%
figure, plot(time,position(5,:),copy.time,copy.position(5,:)), ylabel('Angular Position [rad]'), legend('WECSim','Copy'), xlim([0 copy.time(endInd)])
figure, plot(time,position(5,:),copy.time,copy.velocity(5,:)), ylabel('Angular Velocity [rad/s]'), legend('WECSim','Copy'), xlim([0 copy.time(endInd)])

figure, plot(copy.time,copy.position(1,:),copy.time,5*sin(copy.position(5,:))), ylabel('Surge Position of COG [m]'), legend('from x','From theta')
figure, plot(copy.time,copy.position(3,:),copy.time,5*cos(copy.position(5,:))-5), ylabel('Heave Position of COG [m]'), legend('from z','From theta')

figure, plot(time,forceExcitation(5,:),copy.time,copy.forceExcitation(5,:)), ylabel('Excitation Torque in Pitch [Nm]'), legend('WECSim','Copy'), xlim([0 copy.time(endInd)])
figure, plot(time,forceRadiationDamping(5,:),copy.time(1:101),copy.forceRadiation(5,:)), ylabel('Radiation Torque in Pitch [Nm]'), legend('WECSim','Copy'), xlim([0 copy.time(endInd)])
figure, plot(time,forceRestoring(5,:),copy.time(1:101),copy.forceRestoring(5,:)), ylabel('Restoring Torque in Pitch [Nm]'), legend('WECSim','Copy'), xlim([0 copy.time(endInd)])


%% Functions
function out = WEC_dynamics(states,time,obj,hydro)
% Inputs:
    % states: 2*nDOF x 1
    % time 1x1
    % obj: structure with fields
        % omega:            1 x nOmega      For the excitation force
        % phase:            1 x nOmega      For the excitation force
        % F:             nDOF x nOmega      For the excitation force
        % spectrum:         1 x nOmega      For the excitation force
        % dw:               1 x nOmega      For the excitation force
        % rampTime:         1 x 1           Amount of time for the ramp on the excitation force
        % g:                1 x 1           Gravitational constant
        % cog:           nDOF x 1           Center of gravity of the flap
        % mass:             1 x 1           Structural mass of flap
        % volume:           1 x 1           Submerged volume of the flap
        % velocity:      nDOF x nTime       All velocity measurements (future values are NaN)
        % time              1 x nTime       Time vector which velocity is saved on
        % convolutionDt:    1 x 1           Time step for the convolution integral
        % convolutionLength:1 x 1           Length of time for convolution integral
        % dof               1 x 1           Number of degrees of freedom
        % rad: A structure with fields:
             % Kr:      nDOF x nDOF x nTime
             % time:        1 x nTime
    % hydo: stucture with fields
        % Khs:     nDOF x nDOF x nHeadings      Hydrostatic stiffness matrix

% seperate position and velocity from state vector
position = states(1:obj.dof);
velocity = states((1+obj.dof):2*obj.dof);

%  Excitation force
forceExcitation = GetExcitationForce(obj,time);
    % forceExcitation = interp1(obj.time,obj.forceExcitation2',time)';

% Restoring force
forceRestoring = hydro.Khs(:,:,1) * obj.rho*obj.g * (position-obj.cog) + ...
                      [0;0;obj.g*obj.mass-obj.rho*obj.g*obj.volume;0;0;0];
    % forceRestoring = interp1(obj.time,obj.forceRestoring2',time)';

% Radiation force
forceRadiation = GetRadiationForce(obj,time);
    % forceRadiation = interp1(obj.time,obj.forceRadiation2',time)';

% Total Force
sumOfForces = forceExcitation - forceRadiation - forceRestoring;

% Total Torque
theta = position(5);
sumOfTorque = getNetTorque(sumOfForces,theta);

% Inertial force
mass = [diag(repmat(obj.mass,1,3)), zeros(3,3);zeros(3,3), diag(repmat(obj.inertia,1,3))];
addedMass = hydro.Ainf(1:obj.dof,1:obj.dof)*obj.rho;
I = mass+addedMass;

% Find acceleration - Linerized inertial force = \ddot \theta (I(5,5) + (I(5,1) + I(1,5))*r + I(1,1)*r^2)
r = 5;
thetaddot = sumOfTorque / (I(5,5) + (I(5,1) + I(1,5))*r + I(1,1)*r^2);

% put back into 6 dof
thetadot = velocity(5);
acceleration = [r*thetaddot*cos(theta) - r*thetadot^2*sin(theta);0;-r*thetaddot*sin(theta)-r*thetadot^2*cos(theta);0;thetaddot;0];

% Save to output
out.acceleration = acceleration;
out.forceExcitation = forceExcitation;
out.forceRestoring = forceRestoring;
out.forceRadiation = forceRadiation;
end
    
function forceExcitation = GetExcitationForce(obj,time)
% Inputs: obj is a structure with feilds
    % omega:    1 x nOmega
    % phase:    1 x nOmega
    % F:     nDOF x nOmega
    % spectrum: 1 x nOmega
    % dw:       1 x nOmega
    % rampTime: 1 x 1
% Output: ForceExcitation: 6 x 1 vector

ramp = .5*(1+cos(pi + pi/obj.rampTime*time))*(time<obj.rampTime) + (time>=obj.rampTime);

forceExcitation = ramp* ...
              real(sum( ...
                  ones(obj.dof,1)*exp(1i*(time*obj.omega+obj.phase)) .* ...
                  (obj.F.*sqrt(2*obj.spectrum.*obj.dw))...
              ,2));
end

function forceRadiation = GetRadiationForce(obj,time)
% Inputs: 
    % obj is a structure with feilds:
        % convolutionDt:  1x1 Time step for the convolution integral
        % convolutionLength: 1x1 Length of time for convolution integral
        % dof 1x1 Number of degrees of freedom
        % rad: A structure with fields:
             % Kr: nDOF x nDOF x nTime
            % time: 1 x nTime
    % time: 1x1 Current time that the radiation force is calculated at
% Outputs:
    % forceRadiation: nDOF x 1

convTime = 0:obj.convolutionDt:obj.convolutionLength;
if time > obj.convolutionLength % If the time is greater than our time vector, start the integration obj.convolutionLength seconds before
    startTime = time-obj.convolutionLength; % Start time of convolution integral
    endInd = length(convTime);
else % If the time is less than out time vector, just start at zero and go to the current time
    startTime = 0;
    [~,  endInd] = min(abs(convTime-time)); % ending index corresponding to length of integration
end

% Compute integrand for each value of tau
integrand = NaN(obj.dof,endInd);         % Initialize integrand
for tau_ind = 1:endInd              % Loop over values of tau
    tau = startTime+convTime(tau_ind);
    v = interp1(obj.time,obj.velocity',tau)';% Interpolate velocity since it is on a different time mesh
    Kr = interpolateMatrix(obj.rad,time-tau)*obj.rho;    % times by rho to dimensionalize
    integrand(:,tau_ind) = Kr(1:obj.dof,1:obj.dof)*v;           % Store the integrand for each value of tau
end
% integrate over tau
forceRadiation = trapz(convTime(1:endInd),integrand,2); % perform the integration
end


function torque = getNetTorque(forces,theta)
r = 5;
    torque = forces(5,:)  + r*forces(1,:).*cos(theta) - r*forces(3,:).*sin(theta);
end


function plotForces(time,WECSimForce,hydroForce)
figure,
    subplot(131), plot(time,WECSimForce(1,:),time,hydroForce(1,:)), title('Surge'), grid, ylabel('Force or Torque [N or Nm]')
    subplot(132), plot(time,WECSimForce(3,:),time,hydroForce(3,:)), title('Heave'), grid
    subplot(133), plot(time,WECSimForce(5,:),time,hydroForce(5,:)), title('Pitch'), grid, xlabel('Time [s]')
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

function Kr = interpolateMatrix(rad,time)
% Inputs:
    % rad: A structure with the fields
        % Kr:   nDOF x nDOF x nTime     Radiation impulse response function
        % time:    1 x nTime            Time vector corresponding to Kr
    % time      1x1:            Time at which the matrix Kr is desired
% Outputs:
    % Kr: nDOF x nDOF       Radiation damping convolution kernal at the desires time

n = size(rad.Kr,1);
m = size(rad.Kr,2);
Kr = NaN(n,m);
    for i = 1:n
        for j = 1:m
            Kr(i,j) = interp1(rad.time,squeeze(rad.Kr(i,j,:)),time);
        end
    end
end

