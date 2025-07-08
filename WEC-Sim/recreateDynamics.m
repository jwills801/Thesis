% The goal of this code to recreate the dynamics in WECSim with a State Space time
% domain formulation. Step one is to ensure we can recreate the forces and
% how they combine for F=Ma, which is not as simple in real time because of
% added mass, but should be straight forward to check in post processing
% since acceleration is already known

%%
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
% M = body(1).hydroForce.hf1.storage.fAddedMass*2+[eye(3,3)*body(1).mass, zeros(3,3); zeros(3,3), diag(body(1).inertia)];

massTimesAcceleration = M*acceleration;

% realforceTotal = forceTotal;
realforceTotal = +forceTotal-forcePTO-forceConstraint;
% realforceTotal = -forceAddedMass;


i = 5; % look at just pitch for simplicity
figure, plot(time,massTimesAcceleration(i,:),time,realforceTotal(i,:)),legend('m*a','sum of Forces'), grid

tmp = massTimesAcceleration - realforceTotal;
figure, plot(time,tmp(i,:)), ylabel('m*a-real total force'), grid

figure, plot(time,forceAddedMass(i,:),time,forceExcitation(i,:),time,forceLinearDamping(i,:),time,forceMorisonAndViscous(i,:),time,forceRadiationDamping(i,:),time,forceRestoring(i,:),time,forcePTO(i,:),time,forceConstraint(i,:))
    legend('Added Mass','Excitation','Linear Damping','Morison and Viscous','Radiation Damping','Restoring','PTO','Constraint')
    ylabel('Torque [Nm]'), xlabel('Time [s]'), grid

return
% Below: Check individual forces

%% Check added mass: % Whats with this 1e3??
forceAddedMassCheck = 1e3*body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,1:6)*acceleration;
%forceAddedMassCheck = body(1).hydroForce.hf1.storage.fAddedMass*acceleration;
figure, plot(time,forceAddedMass(5,:),time,forceAddedMassCheck(5,:))
xlabel('Time'), ylabel('Added Mass [Nm]'), legend('WEC-Sim','My Check'), grid

%% Check Linear Damping
figure, plot(time,forceLinearDamping)
xlabel('Time'), ylabel('Linear Damping [N] or [Nm]')

%% Check Morison and Viscous
figure, plot(time,forceMorisonAndViscous)
xlabel('Time'), ylabel('Morison and Viscious Damping [N] or [Nm]')

%% Check Radiation Damping
% resample impulse function to be on the same time scale as velocity (and
% pad extra values with zeros
% Again, what's with the 1e3 factor?
Kr_coarse = 1e3*cat(3,body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.K(:,1:6,:),zeros(6,6,1));
time_coarse = [body(1).hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.t;time(end)];
Kr = NaN(6,6,length(time));
for i = 1:6
    for j = 1:6
        Kr(i,j,:) = interp1(time_coarse, squeeze(Kr_coarse(i,j,:)), time);
    end
end

% Calculate radiation damping at every time
forceRadiationDampingCheck = NaN(size(forceRadiationDamping));
for timeInd = 1:length(time)

tmp = 0;
    for tauInd = 1:timeInd-1
        tmp = tmp + Kr(:,:,timeInd-tauInd)*velocity(:,tauInd)*time(2);
    end
forceRadiationDampingCheck(:,timeInd) = tmp;
end

figure, plot(time,forceRadiationDamping), ylabel('Radiation Damping')
figure, plot(time,forceRadiationDampingCheck), ylabel('Radiation Damping Check')

figure, plot(time,forceRadiationDamping-forceRadiationDampingCheck), ylabel('Damping Difference')

%% Check Restoring Force
forceRestoringCheck = 1e4*body(1).hydroData.hydro_coeffs.linear_restoring_stiffness*position + [simu.rho*simu.gravity*body(1).volume*(body(1).centerBuoyancy-body(1).centerGravity); zeros(3,1)];
figure, plot(time,forceRestoring),ylabel('Restoring'), legend
figure, plot(time,forceRestoringCheck), ylabel('Check Resoring')
figure, plot(time,forceRestoring - forceRestoringCheck), ylabel('Resoring Difference')

%%
% Assemble design matrix X (N timesteps × k forces)
X = [forceExcitation(i,:)', forceAddedMass(i,:)', forceLinearDamping(i,:)', forceMorisonAndViscous(i,:)', forceRadiationDamping(i,:)', forceRestoring(i,:)', forcePTO(i,:)', forceConstraint(i,:)'];
X = [forceRadiationDamping(i,:)', forceRestoring(i,:)', forceConstraint(i,:)'];

% Target vector Y (M·a)
M = [eye(3,3)*body(1).mass, zeros(3,3); zeros(3,3), diag(body(1).inertia)] + body(1).hydroData.hydro_coeffs.added_mass.inf_freq(:,1:6);
massTimesAcceleration = M*acceleration;
Y = massTimesAcceleration(i,:)';

% Solve least squares
W = X \ Y

figure, plot(time,Y,time,X*W)







