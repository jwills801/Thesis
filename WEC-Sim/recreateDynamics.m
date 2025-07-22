% The goal of this code to recreate the dynamics in WECSim with a State Space time
% domain formulation. Step one is to ensure we can recreate the forces and
% how they combine for F=Ma, which is not as simple in real time because of
% added mass, but should be straight forward to check in post processing
% since acceleration is already known

%%
r = 5; % Distance from hinge to center of gravity
time = output.bodies(1).time;

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
figure, plot(time,forceTotal-forceTotalCheck)

%% Check Sum of Forces Equals Mass Time Acceleration
forceConstraint = output.constraints(1).forceConstraint';
forcePTO = output.ptos(1).forceTotal';
% figure, plot(time,forcePTO), ylabel('PTO force')

M = [eye(3,3)*body(1).mass, zeros(3,3); zeros(3,3), diag(body(1).inertia)];
A = body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,1:6)*1e3;
    % Y = sum(diag(body(1).hydroData.hydro_coeffs.added_mass.inf_freq(1:3,1:3)*1e3));
    % dM = [ body(1).adjMassFactor*Y*eye(3), zeros(3,3); zeros(3,3), body(1).hydroData.hydro_coeffs.added_mass.inf_freq(4:6,4:6)*1e3];
    % M_adjusted = M+dM;
    % A_adjusted = A-dM; % Same as body(1).hydroForce.hf1.storage.hydroForce_fAddedMass

massTimesAcceleration = (M)*acceleration;
figure, plot(time,massTimesAcceleration), grid, legend
%%
sumOfForces = forceTotal + forcePTO;

netTorque = sumOfForces(5,:)  + r*sumOfForces(1,:).*cos(position(5,:)) - r*sumOfForces(3,:).*sin(position(5,:));
AccelerationTorque = massTimesAcceleration(5,:)  + r*massTimesAcceleration(1,:).*cos(position(5,:)) - r*massTimesAcceleration(3,:).*sin(position(5,:));
figure, plot(time,massTimesAcceleration(5,:),time,netTorque),legend('m*a','sum of Forces'), grid, xlabel('Time [s]'), ylabel('Torque [Nm]')
figure, plot(time,AccelerationTorque,time,netTorque),legend('m*a','sum of Forces'), grid, xlabel('Time [s]'), ylabel('Torque [Nm]')

figure, plot(time,massTimesAcceleration(5,:) - netTorque), ylabel('m*a-real total force'), grid

i = 5;
figure, plot(time,forceAddedMass(i,:),time,forceExcitation(i,:),time,forceLinearDamping(i,:),time,forceMorisonAndViscous(i,:),time,forceRadiationDamping(i,:),time,forceRestoring(i,:),time,forcePTO(i,:),time,forceConstraint(i,:))
    legend('Added Mass','Excitation','Linear Damping','Morison and Viscous','Radiation Damping','Restoring','PTO','Constraint')
    ylabel('Torque [Nm]'), xlabel('Time [s]'), grid

%% Check added mass:
forceAddedMassCheck = 1e3*body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,1:6)*acceleration;
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
forceRestoringCheck = 1e4*body(1).hydroData.hydro_coeffs.linear_restoring_stiffness*position + [simu.rho*simu.gravity*body(1).volume*(body(1).centerBuoyancy-body(1).centerGravity); zeros(3,1)];
figure, plot(time,forceRestoring),ylabel('Restoring'), legend
figure, plot(time,forceRestoringCheck), ylabel('Check Resoring')
figure, plot(time,forceRestoring - forceRestoringCheck), ylabel('Resoring Difference')

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









