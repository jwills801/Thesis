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