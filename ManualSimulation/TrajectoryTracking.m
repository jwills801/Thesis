
tic
params = getParameters;
cost = dynamics(params);
toc

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
PressureRails = [0 10]*1e6;
rodArea = 0.15^2*pi; % m^2: Radius squared times pi
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
B = [1/(I+Iinf);0;0;0]; C = [1,0,0,0;0,1,0,0];
sys = ss(A,B,C,0);

% Output
params = struct;
params.plotFigures = 1;
params.finalTime = 100;
params.sys = sys;
params.controlInputSet = ptoTorqueOptions(:);

% Excitation
params.period = 5; % s
params.rampTime = 20; % s
end

function mechPower = dynamics(params)
% Time
dt = 1e-3;
t = (0:dt:params.finalTime)';

% Excitation
TexcMag = 1.5e6;
[~,rampStartInd] = min(abs(t-params.rampTime));
ramp = 0.5*(1+cos(pi+pi*t/params.rampTime)).*(t<params.rampTime) + (t>=params.rampTime);
Texc = TexcMag*sin(t*2*pi/params.period).*ramp;

% PI control
H = freqresp(params.sys,2*pi/params.period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');

% Trajectory
thetaDotOpt = Texc/2/real(1/H(1));
thetaOpt = cumtrapz(t,thetaDotOpt);
thetaDDotOpt = [0;diff(thetaDotOpt)./diff(t)];
lambda = 1;
k = 1e-5;
    
% Simulate
states = NaN(length(params.sys.A),length(t)); 
% u = NaN(1,length(t));
states(:,1) = zeros(length(params.sys.A),1);
for timeInd = 1:length(t)-1
    % Control
    thetaDot(timeInd) = states(1,timeInd);
    theta(timeInd) = states(2,timeInd);


    % PI control
        % u(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));
    
    % Sliding Mode Contrl
    thetaError = theta(timeInd)-thetaOpt(timeInd);
    thetaDotError = thetaDot(timeInd)-thetaDotOpt(timeInd);
        % Estimate ThetaDDot without u
    statesDotHat = params.sys.A*states(:,timeInd) + params.sys.B*Texc(timeInd);
    thetaDDotError = statesDotHat(1) - thetaDDotOpt(timeInd);
    
    s = thetaDotError+lambda*thetaError;
        % Set sdot to -k*sign(s)
    u(timeInd) = (-thetaDDotError - lambda*thetaDotError - k*sign(s))/params.sys.B(1);
    
    % Select closest force 
    % [~,controlInd] = min(abs(params.controlInputSet-u(timeInd)));
    % u(timeInd) = params.controlInputSet(controlInd);

    
    if s > 0
        % Pick the largest control option still smaller than u
        smallerOptions = params.controlInputSet(find(params.controlInputSet<=u(timeInd)));
        if ~isempty(smallerOptions)
            u(timeInd) = max(smallerOptions);
        else
            u(timeInd) = min(params.controlInputSet);
        end
    else
        % Pick the smallest option still larger than u
        largerOptions = params.controlInputSet(find(params.controlInputSet>=u(timeInd)));
        if ~isempty(largerOptions)
            u(timeInd) = min(largerOptions);
        else
            u(timeInd) = max(params.controlInputSet);
        end
 
    end


    % Forward Euler Step
    states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(u(timeInd)+Texc(timeInd)))*dt;
end
% For the last time step
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
u(timeInd+1) = -1*(Kp*thetaDot(timeInd+1) + Ki*theta(timeInd+1));
    
% Absorbed Power
mechPower = -u.*thetaDot;
mechEnergy = cumtrapz(t,mechPower);
avePow = trapz(t(rampStartInd:end),mechPower(rampStartInd:end))/(params.finalTime-params.rampTime)

% Theoretical max
avePowOpt = TexcMag^2/8/real(1/H(1))

disp([num2str((avePowOpt-avePow)/avePowOpt*100),'% from optimal'])

if params.plotFigures
    figure
    subplot(221), plot(t,u), xlabel('Time [s]'), ylabel('Control Input [Nm]')
    subplot(222), yyaxis left, plot(t,thetaDot), ylabel('Angular Velocity [rad/s]')
               yyaxis right,plot(t,Texc), xlabel('Time [s]'), ylabel('Excitaiton Torque [Nm]')
    subplot(223), plot(t,thetaDot,'k',t,thetaDotOpt,'k--'), xlabel('Time [s]'), ylabel('Angular Velocity [W]'), legend('Actual','Optimal')
    subplot(224), plot(t,mechPower), xlabel('Time [s]'), ylabel('Absorbed Energy [J]')
end

end
