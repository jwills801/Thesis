% clear, close all
tic
%%
Params = getParameters;

Control = dynamics(Params);
%%
SwitchMap = MakeSwitchingLossMap(Params);
%%
Loss = UseSwitchingLossMap(Params,Control,SwitchMap);

toc

function out = UseSwitchingLossMap(Params,Control,SwitchMap)
    % Calculate volume in each side of the cylinder 
    d = 7;
    r_cyl = 3;
    L_equilib = sqrt(d^2+r_cyl^2 + 2*d*r_cyl*sin(0));
    L_extra = L_equilib - Params.stroke/2; % All extra length that isnt cap volume length
    L = sqrt(d^2+r_cyl^2 + 2*d*r_cyl*sin(Control.theta));
    xCap = L-L_extra;
    xRod = Params.stroke - xCap;
    capVol = xCap*Params.capArea;
    rodVol = xRod*Params.rodArea;

    % Calculate velocity of the cylinder
    dLdt = d*r_cyl*cos(Control.theta).*Control.thetaDot./L;
    capVelA = dLdt*Params.capArea;
    rodVelA = -dLdt*Params.rodArea;

% Find switching events
    % These denote the start of the switch
    % The start and end of the simulation are also counted as events
eventInds = [1, find(diff(Control.u)~=0), length(Control.t)-1];
eventTimes = Control.t(eventInds);

% Previous and new control values
switchFrom = Control.u(eventInds);
switchTo = Control.u(eventInds+1);

mapDT = .1;

% Initilize the energy loss vector
    % This denotes the energy lost between event times
    % Thus thereis one less loss entry than there are events
loss = NaN(length(eventInds)-1,1);
for k = 1:length(eventInds)-1
    % penalize switch
        % Which pressure rail did we switch from
        % ptoTorqueOptions is a matrix
            % Each row is a different cap side option
            % Each col is a different rod side option
    [capSwitchFromInd,rodSwitchFromInd] = find(switchFrom(k)==Params.controlInputSet);
    
        % Which pressure rails are we switching to
    [capSwitchToInd,rodSwitchToInd] = find(switchTo(k)==Params.controlInputSet);

    capSwitchFrom = Params.PR(capSwitchFromInd);
    capSwitchTo = Params.PR(capSwitchToInd);
    capSwitchVelA = capVelA(eventInds(k));
    capSwitchVol = capVol(eventInds(k));
    capSwitchingLoss = interpn(SwitchMap.PR,SwitchMap.PR,SwitchMap.velA_vals,SwitchMap.vol_vals, SwitchMap.Eloss,...
        capSwitchFrom,capSwitchTo,capSwitchVelA,capSwitchVol);

    rodSwitchFrom = Params.PR(rodSwitchFromInd); 
    rodSwitchTo = Params.PR(rodSwitchToInd);
    rodSwitchVelA = rodVelA(eventInds(k));
    rodSwitchVol = rodVol(eventInds(k));
    rodSwitchingLoss = interpn(SwitchMap.PR,SwitchMap.PR,SwitchMap.velA_vals,SwitchMap.vol_vals, SwitchMap.Eloss,...
        rodSwitchFrom,rodSwitchTo,rodSwitchVelA,rodSwitchVol);


    % Penalize time between switches
        % What is the flow over this time period 
    ind1 = eventInds(k) + mapDT/(Control.t(2)-Control.t(1));
    ind2 = eventInds(k+1);
    % Find open valve loss from ind1 to ind2
    capPowerLoss = (abs(capVelA(ind1:ind2))).^3/(SwitchMap.valveConstant^2);
    rodPowerLoss = (abs(rodVelA(ind1:ind2))).^3/(SwitchMap.valveConstant^2);
    steadyLoss = trapz(Control.t(ind1:ind2),capPowerLoss+rodPowerLoss);

    loss(k) = capSwitchingLoss + rodSwitchingLoss + steadyLoss;
end
TotalValveLoss = sum(loss);

% Plot energy loss per event
figure, yyaxis right
plot(eventTimes(1:end-1),loss,'*')
yyaxis left
plot(Control.t,Control.u)

% Plot cummulative valve loss over time
figure, plot(eventTimes(1:end-1),cumsum(loss))

% Output results
out.TotalValveLoss = TotalValveLoss;
end

function params = getParameters(~)

% Stiffness
rho = 1024;
Vol = 297;
g = 9.81;
r_cob = 4;
r_cog = 5;
m = 127000;
Khs = rho*Vol*g*r_cob - m*g*r_cog;

% Inertia (and added inertia)
I = 5.025e6;

% Radiation
Iinf = 1.734e7;
x_timeDomain = [1.9754, 1.1345, 7.6921];

% Set Valued Control Inputs
PressureRails = [0 30]*1e6;
stroke = 5;
rodArea = 0.2^2*pi; % m^2: Radius squared times pi
capArea = 1.5*rodArea; % m^2: Area ratio times rod Area
r= 1.18; % Moment arm from force to torque
capTorqueOptions = PressureRails*capArea*r;
rodTorqueOptions = PressureRails*rodArea*r;
ptoTorqueOptions = capTorqueOptions'-rodTorqueOptions;


% State space
A = [0,            -Khs/(I+Iinf),0,-1/(I+Iinf);...
    1,               0, 0,        0;...
    0,               0, 0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0, 1, -x_timeDomain(2)];
B = [1/(I+Iinf);0;0;0]; C = [1,0,0,0]; C = eye(4);
sys = ss(A,B,C,0);

% Output
params = struct;
params.plotFigures = 1;
params.finalTime = 50;
params.sys = sys;

params.PR = PressureRails;
params.controlInputSet = ptoTorqueOptions;
params.capArea = capArea;
params.rodArea = rodArea;
params.stroke = stroke;

% Excitation
params.period = 10; % s
params.rampTime = 1e-3; % s
end

function out = dynamics(params)
% Time
dt = 1e-2;
t = (0:dt:params.finalTime)';

% Control Horizon
controlHorizon = 2e-1;

% How many simulation time steps are there in a control horizon
HorizonSampling = round(controlHorizon/dt);  % Update every N simulation steps 

% Excitation
% Regular Waves
% TexcMag = 2e6;
% [~,rampStartInd] = min(abs(t-params.rampTime));
% ramp = 0.5*(1+cos(pi+pi*t/params.rampTime)).*(t<params.rampTime) + (t>=params.rampTime);
% Texc = TexcMag*sin(t*2*pi/params.period).*ramp;

% Irregular Waves
Tp = params.period;
Hs = 1;
rampTime = params.rampTime;
finalTime = params.finalTime;
irregWaves = IrregularWaves(Tp,Hs,rampTime,finalTime,dt);
Texc = irregWaves.Texc;
[~,rampStartInd] = min(abs(t-params.rampTime));

Z = freqresp(params.sys,irregWaves.waveInfo.omega);
Rr = real(1./squeeze(Z(1,1,:)))';
avePowOpt = trapz(irregWaves.waveInfo.omega,abs(irregWaves.F).^2.*irregWaves.waveInfo.spectrum*2/8./Rr)

% PI control
H = freqresp(params.sys,2*pi/params.period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');

% Trajectory
ramp = .5*(1+cos(pi + pi/rampTime*t)).*(t<rampTime) + (t>=rampTime);
rng(1) % Sets random seed for repeatability
phase = 2*pi*rand(size(irregWaves.waveInfo.omega));
dw = [0,diff(irregWaves.waveInfo.omega)];

thetaDotOpt = ramp.* real(sum( exp(1i*(t*irregWaves.waveInfo.omega+phase)) .* ...
    (irregWaves.F/2./Rr.*sqrt(2*irregWaves.waveInfo.spectrum.*dw)) ,2));

% thetaDotOpt = Texc/2/real(1/H(1));
thetaOpt = cumtrapz(t,thetaDotOpt);
thetaDDotOpt = [0;diff(thetaDotOpt)./diff(t)];
lambda = 1;
k = 1e-5;
phi = 3e-2; % band around sliding surface
    
% Simulate
states = NaN(length(params.sys.A),length(t)); 
u = NaN(1,length(t)); u(1) = 0;
states(:,1) = zeros(length(params.sys.A),1);

waitbarObj = waitbar(0,'Simulating WEC Dynamics');
for timeInd = 1:length(t)-1
    waitbar(timeInd/length(t),waitbarObj);

    thetaDot(timeInd) = states(1,timeInd);
    theta(timeInd) = states(2,timeInd);

    % Sliding Mode Contrl
    thetaError = theta(timeInd)-thetaOpt(timeInd);
    thetaDotError = thetaDot(timeInd)-thetaDotOpt(timeInd);

    % Estimate ThetaDDot without u
    statesDotHat = params.sys.A*states(:,timeInd) + params.sys.B*Texc(timeInd);
    thetaDDotError = statesDotHat(1) - thetaDDotOpt(timeInd);

    s(timeInd) = thetaDotError+lambda*thetaError;

    % if s is large, recalulate control input
    if abs(s(timeInd)) > phi || timeInd == 1
        % Horizon indices
        hInds = (timeInd:min((timeInd+HorizonSampling-1),length(t)));

        xStar = [thetaOpt(hInds),thetaDotOpt(hInds),zeros(length(hInds),2)];

        lastSwitchInd = max([find(diff(u(1:max(timeInd-1,1)))~=0,1,'last'),1]);
        timeSinceLastSwitch = t(timeInd)-t(lastSwitchInd);
        previousInput = find(u(lastSwitchInd+1)==params.controlInputSet(:));
        J = 1*exp(-10*timeSinceLastSwitch) * ones(size(params.controlInputSet(:)));
        J(previousInput) = 0;

        % Simulate future times for each control input
        for uInd = 1:length(params.controlInputSet(:))

            % Control inputs:
            uOption = params.controlInputSet(uInd);

            xHat = lsim(params.sys,uOption+Texc(hInds),t(hInds),states(:,timeInd));

            S = (xHat - xStar)*[1; lambda; 0; 0];

            % J(uInd) = J(uInd) + sum(S.^2);
            J(uInd) = J(uInd) + abs(S(end));
        end
        [~,uOptInd] = min(J);
   
        u(timeInd) = params.controlInputSet(uOptInd);


    else % If s is small, use previous control input
        u(timeInd) = u(timeInd-1);
    end % if statement on size of s

    % PI control
    % u(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));

    % Forward Euler Step
    % states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(u(timeInd)+Texc(timeInd)))*dt;
    
    % Forward with lsim
    stepInds = [timeInd,timeInd+1];
    statesMany = lsim(params.sys,u(timeInd)+Texc(stepInds),t(stepInds),states(:,timeInd));
    states(:,timeInd + 1) = statesMany(end,:);
end
close(waitbarObj)

% For the last time step
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
u(timeInd+1) = u(timeInd);
    
% Absorbed Power
mechPower = -u.*thetaDot;
figure, plot(t,cumtrapz(mechPower))
mechEnergy = cumtrapz(t,mechPower);
avePow = trapz(t(rampStartInd:end),mechPower(rampStartInd:end))/(params.finalTime-params.rampTime)

% Theoretical max
% avePowOpt = 0; % TexcMag^2/8/real(1/H(1))

% Count switching events
switchInds = find(diff(u) ~= 0) + 1;  % +1 to get the index of the new value
switchRate = length(switchInds)/finalTime

%% Output
out = struct();
out.u = u;
out.t = t;
out.theta = theta;
out.thetaDot = thetaDot;
out.optimalAvergePower = avePowOpt;
out.averagePower = avePow;
out.mechanicalPower = mechPower;
out.switchRate = switchRate;

%% Display results

disp([num2str((avePowOpt-avePow)/avePowOpt*100),'% from optimal'])

if params.plotFigures
    figure
    subplot(221), plot(t,u), xlabel('Time [s]'), ylabel('Control Input [Nm]'), grid
    subplot(222), yyaxis left, plot(t,thetaDot), ylabel('Angular Velocity [rad/s]'), grid
               yyaxis right,plot(t,Texc), xlabel('Time [s]'), ylabel('Excitaiton Torque [Nm]')
    subplot(223), plot(t,thetaDot,'k',t,thetaDotOpt,'k--'), xlabel('Time [s]'), ylabel('Angular Velocity [W]'), legend('Actual','Optimal'), grid
    subplot(224), plot(t,[0,s]), xlabel('Time [s]'), ylabel('Sliding Surface'), grid
end

end
