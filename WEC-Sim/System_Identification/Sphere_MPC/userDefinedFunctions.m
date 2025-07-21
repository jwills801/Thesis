%Example of user input MATLAB file for post processing
close all

%Plot waves
waves.plotElevation(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end

%Plot heave response for body 1
output.plotResponse(1,3);

%Plot heave forces for body 1
output.plotForces(1,3);

controllersOutput = controller1_out;
signals = {'force','power'};
for ii = 1:length(controllersOutput)
    for jj = 1:length(signals)
        controllersOutput(ii).(signals{jj}) =  controllersOutput(ii).signals.values(:,(jj-1)*6+1:(jj-1)*6+6);
    end
end

% Plot Controller Power
figure()
plot(controllersOutput.time,controllersOutput.power(:,3)/1000)
title('Controller Power')
ylabel('Power (kW)')
xlabel('Time (s)')

% Calculate average controller power over last 10 wave periods
endInd = length(controllersOutput.power(:,3));
startTime = controllersOutput.time(end) - 10*waves.period; % select last 10 periods
[~,startInd] = min(abs(controllersOutput.time(:) - startTime));
disp('Controller Power (kW):')
disp(mean(controllersOutput.power(startInd:endInd,3))/1000)
disp("Peak Controller Power (kW):")
disp(max(abs(controllersOutput.power(startInd:endInd,3)))/1000)

% Plot heave position
figure()
plot(output.bodies.time,output.bodies.position(:,3)-body.centerGravity(3))
title('Heave Position')
yline(controller.modelPredictiveControl.maxPos, '--')
yline(-controller.modelPredictiveControl.maxPos, '--')
ylabel('Position (m)')
xlabel('Time (s)')
ylim([-(controller.modelPredictiveControl.maxPos + .25), controller.modelPredictiveControl.maxPos + .25])
disp("Max Heave Position (m):")
disp(max(abs(output.bodies.position(:,3)-body.centerGravity(3))))

% Plot heave velocity
figure()
plot(output.bodies.time,output.bodies.velocity(:,3))
title('Heave Velocity')
yline(controller.modelPredictiveControl.maxVel, '--')
yline(-controller.modelPredictiveControl.maxVel, '--')
ylabel('Velocity (m/s)')
xlabel('Time (s)')
ylim([-(controller.modelPredictiveControl.maxVel + .25), controller.modelPredictiveControl.maxVel + .25])
disp("Max Heave Velocity (m/s):")
disp(max(abs(output.bodies.velocity(:,3))))

% Plot PTO Force
figure()
plot(controllersOutput.time,controllersOutput.force(:,3)/1000)
title('PTO Force')
yline(controller.modelPredictiveControl.maxPTOForce/1000, '--')
yline(-controller.modelPredictiveControl.maxPTOForce/1000, '--')
ylabel('Force (kN)')
xlabel('Time (s)')
ylim([-(controller.modelPredictiveControl.maxPTOForce/1000 + 250), controller.modelPredictiveControl.maxPTOForce/1000 + 250])
disp("Max PTO Force (kN):")
disp(max(abs(controllersOutput.force(:,3)))/1000)

% Plot change in PTO Force
figure()
plot(controllersOutput.time,gradient(controllersOutput.force(:,3))*(1/simu.dt)/1000)
title('PTO Force Change')
yline(controller.modelPredictiveControl.maxPTOForceChange/1000, '--')
yline(-controller.modelPredictiveControl.maxPTOForceChange/1000, '--')
ylabel('Force Change (kN/s)')
xlabel('Time (s)')
ylim([-(controller.modelPredictiveControl.maxPTOForceChange/1000 + 250), controller.modelPredictiveControl.maxPTOForceChange/1000 + 250])
disp("Max PTO Force Change (kN/s):")
disp(max(abs(gradient(controllersOutput.force(:,3))*(1/simu.dt)))/1000)

%% Simulate the system using the system matrices to make sure it matches WEC-SIM results
dt = time(2);
dt_d = .5;
sys_d = c2d(controller(1).plant.sys_c,dt_d,'zoh');
time = output.bodies.time;
zb = output.bodies.position(:,3)-body.centerGravity(3);

Fe = output.bodies(1).forceExcitation(:,3);
Fe_d = Fe(1:dt_d/dt:end);
time_d = time(1:dt_d/dt:end);

F_PTO = output.ptos(1).forceTotal(:,3);
dF_PTO = [0;diff(F_PTO)]/time(2);
dF_PTO_d = dF_PTO(1:dt_d/dt:end);

X = NaN(length(time_d),length(sys_d.A));
X(1,:) = zeros(1,length(sys_d.A));
for i = 2:length(time_d)
    X(i,:) = sys_d.A*X(i-1,:)' + sys_d.B*[dF_PTO_d(i-1);Fe_d(i-1)];
end
dzb_check = X(:,1);
zb_check = X(:,2);
F_PTO_check = X(:,9);

figure
subplot(311), plot(time_d,dzb_check,time,output.bodies.velocity(:,3))
subplot(312), plot(time_d,zb_check,time,zb)
subplot(313), plot(time_d,F_PTO_check,time,F_PTO)

%% Simulate the continuous time system
sys_c = controller(1).plant.sys_c;
X = NaN(length(time),length(sys_c.A));
X(1,:) = zeros(1,length(sys_c.A));
dX = X;
for i = 2:length(time)
    dX(i,:) = sys_c.A*X(i-1,:)' + sys_c.B*[dF_PTO(i-1);Fe(i-1)];
    X(i,:) = X(i-1,:)' + dX(i,:)'*dt;
end
dzb_check = X(:,1);
zb_check = X(:,2);
F_PTO_check = X(:,9);

figure
subplot(311), plot(time,dzb_check,time,output.bodies.velocity(:,3))
subplot(312), plot(time,zb_check,time,zb)
subplot(313), plot(time,F_PTO_check,time,F_PTO)

%%
acc = [0;diff(dzb_check)]/dt;
figure, plot(time,acc,time,dX(:,1))

%% Sum of Forces Equals Mass Times Acceleration
Ma = (controller(1).bemData.m+controller(1).bemData.aInf)*dX(:,1);
sumOfForces = output.bodies(1).forceTotal(:,3);
sumOfForces = output.bodies(1).forceTotal(:,3) + F_PTO + output.bodies(1).forceAddedMass(:,3);
figure, plot(time,Ma,time,sumOfForces)



