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
% figure, plot(hydro.w,spectrum), ylabel('Spectrum [m^2 s/rad]') % plotSpectrum(waves)

hydro.forceExcitation = NaN(size(hydro.forceRestoring));
for i = 1:sum(hydro.dof)
    ex_re = interp1(hydro.w,squeeze(hydro.ex_re(i,1,:)),omega);
    ex_im = interp1(hydro.w,squeeze(hydro.ex_im(i,1,:)),omega);
    F = squeeze((ex_re'+ex_im'*1i)) * rho*g;
    hydro.forceExcitation(i,:) = real(exp(1i*(time*omega+phase)) * (F.*sqrt(2*spectrum.*dw)));
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
dof = 1;
figure, plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:)))

dof = 3;
figure, plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:)))

dof = 5;
figure, plot(hydro.ra_t,squeeze(hydro.ra_K(dof,dof,:)),rad.time,squeeze(rad.Kr(dof,dof,:)),rad.time,squeeze(rad.Kr2(dof,dof,:)))

dof = 6; figure, plot(hydro.w,squeeze(hydro.B(dof,dof,:)),body(1).hydroData.simulation_parameters.w,squeeze(body(1).hydroData.hydro_coeffs.radiation_damping.all(dof,dof,:)))

%%    % Convolution Integral
tic
rad.convTime = 0:.1:30;
rad.convolutionLength = 10; % It takes about 10 seconds for Kr terms to decay plot(hydro.w,reshape(hydro.B,[12*12,500]))
F = NaN(12,length(time));   % Initilize force
for time_ind = 1:301%length(time)         % Loop over time
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
% figure, plot(time,F)
figure, plot(time,F(5,:),time,forceRadiationDamping(5,:))
toc
return

%% Check Sum of Forces Equals Mass Time Acceleration
forceConstraint = output.constraints(1).forceConstraint';
forcePTO = output.ptos(1).forceTotal';
% figure, plot(time,forcePTO), ylabel('PTO force')

M = [eye(3,3)*body(1).mass, zeros(3,3); zeros(3,3), diag(body(1).inertia)];
A = body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,1:6)*simu.rho;

massTimesAcceleration = (M+A)*acceleration;
accelerationTorque = massTimesAcceleration(5,:)  + r*massTimesAcceleration(1,:).*cos(position(5,:)) - r*massTimesAcceleration(3,:).*sin(position(5,:));
    % Linearized version - found in symbolicLinearizaton.m: theta_ddot*(I_55 + I_51*r + r*(I_15 + I_11*r))
    I = M+A;
    accelerationTorqueLinearized = ( I(5,5) + r*(I(5,1) + I(1,5) + r*I(1,1)) )*acceleration(5,:);
figure, plot(time,accelerationTorque,time,accelerationTorqueLinearized)
%%
sumOfForces = forceTotal + forceAddedMass;
netTorque = sumOfForces(5,:)  + r*sumOfForces(1,:).*cos(position(5,:)) - r*sumOfForces(3,:).*sin(position(5,:));

figure, plot(time,accelerationTorque,time,netTorque),legend('m*a','sum of Forces'), grid, xlabel('Time [s]'), ylabel('Torque [Nm]')

figure, plot(time,accelerationTorque - netTorque), ylabel('m*a-real total force'), grid

i = 5;
figure, plot(time,forceAddedMass(i,:),time,forceExcitation(i,:),time,forceLinearDamping(i,:),time,forceMorisonAndViscous(i,:),time,forceRadiationDamping(i,:),time,forceRestoring(i,:),time,forcePTO(i,:),time,forceConstraint(i,:))
    legend('Added Mass','Excitation','Linear Damping','Morison and Viscous','Radiation Damping','Restoring','PTO','Constraint')
    ylabel('Torque [Nm]'), xlabel('Time [s]'), grid


theta = position(5,:);
figure, plot(time,getNetTorque(forceAddedMass,theta),time,getNetTorque(forceExcitation,theta),time,getNetTorque(forceLinearDamping,theta),time,getNetTorque(forceMorisonAndViscous,theta),time,getNetTorque(forceRadiationDamping,theta),time,getNetTorque(forceRestoring,theta),time,getNetTorque(forcePTO,theta),time,getNetTorque(forceConstraint,theta))
    legend('Added Mass','Excitation','Linear Damping','Morison and Viscous','Radiation Damping','Restoring','PTO','Constraint')
    ylabel('Torque [Nm]'), xlabel('Time [s]'), grid

%% Check added mass:
forceAddedMassCheck = A*acceleration;
figure, plot(time,forceAddedMass(5,:),time,forceAddedMassCheck(5,:))
xlabel('Time'), ylabel('Added Mass [Nm]'), legend('WEC-Sim','My Check'), grid

%% Check Linear Damping
figure, plot(time,forceLinearDamping)
xlabel('Time'), ylabel('Linear Damping [N] or [Nm]')

%% Check Morison and Viscous
figure, plot(time,forceMorisonAndViscous)
xlabel('Time'), ylabel('Morison and Viscious Damping [N] or [Nm]')

%% Check Radiation Damping
%% Using the convolution
% % resample impulse function to be on the same time scale as velocity (and
% % pad extra values with zeros
% Kr_coarse = 1e3*cat(3,body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.K(:,1:6,:),zeros(6,6,1));
% time_coarse = [body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.t;time(end)];
% Kr = NaN(6,6,length(time));
% for i = 1:6
%     for j = 1:6
%         Kr(i,j,:) = interp1(time_coarse, squeeze(Kr_coarse(i,j,:)), time);
%     end
% end
% 
% % Calculate radiation damping at every time
% forceRadiationDampingCheck = NaN(size(forceRadiationDamping));
% for timeInd = 1:length(time)
% 
% tmp = 0;
%     for tauInd = 1:timeInd-1
%         tmp = tmp + Kr(:,:,timeInd-tauInd)*([1;1;1;1;1/r;1].*velocity(:,tauInd))*time(2);
%     end
% forceRadiationDampingCheck(:,timeInd) = tmp;
% end
% 
% figure, plot(time,forceRadiationDamping), ylabel('Radiation Damping'), grid, legend
% figure, plot(time,forceRadiationDampingCheck), ylabel('Radiation Damping Check: Convolution'), grid, legend
% figure, plot(time,forceRadiationDamping-forceRadiationDampingCheck), ylabel('Damping Difference'), grid, legend

%% Using State space Realization of radiation damping
sys_rad = ss(body(1).hydroForce.hf1.ssRadf.A,body(1).hydroForce.hf1.ssRadf.B,body(1).hydroForce.hf1.ssRadf.C,body(1).hydroForce.hf1.ssRadf.D);
forceRadiationDampingCheck = lsim(sys_rad,velocity,time);
% figure, plot(time,forceRadiationDamping), ylabel('Radiation Damping'), grid, legend
% figure, plot(time,forceRadiationDampingCheck), ylabel('Radiation Damping Check: SS'), grid, legend
figure, plot(time,forceRadiationDamping-forceRadiationDampingCheck'), ylabel('Damping Difference'), grid, legend


%% Check Restoring Force
forceRestoringCheck = body(1).hydroForce.hf1.linearHydroRestCoef*(position-[body(1).centerGravity;0;0;0]) + ...
    [zeros(2,1);simu.gravity*body(1).mass-simu.rho*simu.gravity*body(1).volume;zeros(3,1)];
figure, plot(time,forceRestoring - forceRestoringCheck), ylabel('Restoring Difference')

K = body(1).hydroForce.hf1.linearHydroRestCoef;
torqueRestoring = K(5,5)*position(5,:) - r*sin(position(5,:)).*(K(3,3)*(r*cos(position(5,:))-r)+simu.gravity*(body(1).mass-simu.rho*body(1).volume));
torqueRestoringLinearized = (K(5,5) - r*simu.gravity*(body(1).mass-simu.rho*body(1).volume))*position(5,:);
figure, plot(time,torqueRestoring,time,torqueRestoringLinearized)


%% Full state space estimation of dynamics
a = squeeze(body(1).hydroData.hydro_coeffs.added_mass.all(5,5,:)).*simu.rho;
aInf = body(1).hydroData.hydro_coeffs.added_mass.inf_freq(5,5)*simu.rho;
m = body(1).mass;
k = body(1).hydroData.hydro_coeffs.linear_restoring_stiffness(5,5)*simu.rho*simu.gravity; 
b = squeeze(body(1).hydroData.hydro_coeffs.radiation_damping.all(5,5,:)).*simu.rho.*body(1).hydroData.simulation_parameters.w';
load('coeff.mat') % load transfer function found in plotFreqDep.m

controller(1).plant.A = zeros(8);
controller(1).plant.A(1,2) = -k/(m+aInf);
controller(1).plant.A(1,5) = 1/(m+aInf);
controller(1).plant.A(2:4,1:3)=eye(3);
controller(1).plant.A(5,:)=[-coeff.KradNum -coeff.KradDen(2:end)];
controller(1).plant.A(6:end,5:end-1) = eye(3);

controller(1).plant.Bv = zeros(8,1);
controller(1).plant.Bv(1,1) = 1/(m+aInf);
controller(1).plant.Bu = zeros(8,1);
controller(1).plant.Bu(1,1) = 1/(m+aInf);

controller(1).plant.C = [1, zeros(1,8-1); ...
    0, 1, zeros(1,8-2)];

% Augmenting for the delta formulation which has Fpto in the states and replaces that w/ dFpto as the control input
controller(1).plant.A = [controller(1).plant.A, controller(1).plant.Bu; zeros(1,8+1)];
controller(1).plant.Bu = [zeros(8,1); 1];
controller(1).plant.Bv = [controller(1).plant.Bv; 0];  
controller(1).plant.C = [controller(1).plant.C zeros(2,1); zeros(1,8), 1]; 

controller(1).plant.Du  = zeros(3,1);   % dFpto has no contripution to output states. [3x1] per each pod since output is Z, dZ, Fpto
controller(1).plant.Dv = zeros(3,1);    % Fe has no contripution to output states

% simulate the system
dF_PTO = [0,diff(forcePTO(5,:))]/time(2);
Fe = forceExcitation(5,:);
sys_c = ss(controller(1).plant.A,[controller(1).plant.Bu controller(1).plant.Bv],controller(1).plant.C,[controller(1).plant.Du controller(1).plant.Dv]);
X = NaN(length(time),length(sys_c.A));
X(1,:) = zeros(1,length(sys_c.A));
dX = X;
for i = 2:length(time)
    dX(i,:) = sys_c.A*X(i-1,:)' + sys_c.B*[dF_PTO(i-1);Fe(i-1)];
    X(i,:) = X(i-1,:)' + dX(i,:)'*time(2);
end
dzb_check = X(:,1);
zb_check = X(:,2);
F_PTO_check = X(:,9);

figure
subplot(311), plot(time,dzb_check,time,velocity(5,:)), ylim([-10 10])
subplot(312), plot(time,zb_check,time,position(5,:))
subplot(313), plot(time,F_PTO_check,time,forcePTO(5,:))



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
g = 9.81; rho = 1000;
frequency = omega'/(2*pi);
dw = [omega(2)-omega(1); diff(omega')];

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
spectrum = fSpectrum/(2*pi);

waveNumber = calcWaveNumber(omega,waterDepth,g,0); %Calculate Wave Number for Larger Number of Frequencies Before Down Sampling in Equal Energy Method
power = sum((1/2)*rho*g*spectrum'.*dw'.*sqrt(g./waveNumber.*tanh(waveNumber.*waterDepth)).*(1 + 2.*waveNumber.*waterDepth./sinh(2.*waveNumber.*waterDepth)));

end

function Kr = interpolateMatrix(rad,time)
n = size(rad.Kr,1);
m = size(rad.Kr,2);
Kr = NaN(n,m);
    for i = 1:n
        for j = 1:m
            Kr(i,j) = interp1(rad.time,squeeze(rad.Kr(i,j,:)),time);
        end
    end
end

