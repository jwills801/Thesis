params = getParameters;
params.plotFigures=1;
params.rounding = 0;
cost = PI(params)

params.rounding = 1;
cost = PI(params)

function params = getParameters(~)
rho = 1024;
Vol = 297;
g = 9.81;
r_cob = 4;
r_cog = 5;
m = 127000;
Khs = rho*Vol*g*r_cob - m*g*r_cog;
I = 5.025e6;
Iinf = 1.734e7;
x_timeDomain = [1.9754, 1.1345, 7.6921];

PR = [0 10]*1e6;
% PR = [0 8 15]*1e6;
% PR = [0 5 10 15]*1e6;
rodArea = 0.15^2*pi; % m^2: Radius squared times pi
capArea = 1.5*rodArea; % m^2: Area ratio times rod Area
r= 1.18; % Moment arm from force to torque
capTorqueOptions = PR*capArea*r;
rodTorqueOptions = PR*rodArea*r;
ptoTorqueOptions = capTorqueOptions'-rodTorqueOptions;
ptoTorqueOptions(4) = NaN; 

% State space
A = [0,            -Khs/(I+Iinf),0,-1/(I+Iinf);...
    1,               0, 0,        0;...
    0,               0, 0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0, 1, -x_timeDomain(2)];
B = [1/(I+Iinf);0;0;0]; C = [1,0,0,0;0,1,0,0];
sys = ss(A,B,C,0);

params = struct;
params.plotFigures = 0;
params.uDT = .1; % optimization step
params.finalTime = 50;
params.ptoTorqueOptions = ptoTorqueOptions;
params.PR = PR;
params.capArea = capArea; params.rodArea = rodArea; params.r = r;

    % load ExcitationTorque.mat % Only gives forces at 0.2s intervals
    % params.T_exc = T_exc(1:length(params.t),2);
params.period = 5; % s
params.rampTime = 20; % s

params.sys = sys; params.r=r; params.capArea=capArea; params.rodArea=rodArea;
end

function cost = PI(params)
dt = 1e-3;
controlStartTime = 25;
t = (0:dt:params.finalTime)';
TorquePTO = NaN(size(t));

H = freqresp(params.sys,2*pi/params.period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');
ramp = 0.5*(1+cos(pi+pi*t/params.rampTime)).*(t<params.rampTime) + (t>=params.rampTime);
Texc = 1.5e6*sin(t*2*pi/params.period).*ramp;
    
states = NaN(length(params.sys.A),length(t)); 
states(:,1) = zeros(length(params.sys.A),1);
for timeInd = 1:length(t)-1
    thetaDot(timeInd) = states(1,timeInd);
    theta(timeInd) = states(2,timeInd);

    TorquePTO(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));
    if params.rounding ==1
        [~,roundInd] = min(abs(TorquePTO(timeInd) - params.ptoTorqueOptions(:)));
        TorquePTO(timeInd) = params.ptoTorqueOptions(roundInd);
    end

    states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(TorquePTO(timeInd)+Texc(timeInd)))*dt;
end
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
TorquePTO(timeInd+1) = -1*(Kp*thetaDot(timeInd+1) + Ki*theta(timeInd+1));
    
mechPower = TorquePTO'.*thetaDot;

cost = sum(mechPower)*dt;

cost = cost/1e6;


if params.plotFigures
    Pmax = (1.5e6)^2/8/Kp;
    Pact = cost/t(end)*1e6;
    Pact_kW = Pact/1e3;
    thetaDotOpt_Theory = Texc/2/Kp;
    thetaOpt_Theory = cumtrapz(t,thetaDotOpt_Theory);
    figure, subplot(211), plot(t,TorquePTO,t,-1*(Kp*thetaDotOpt_Theory + Ki*thetaOpt_Theory)), ylabel('PTO Torque [Nm]'), xlabel('Time [s]'), %legend('Numerical Optimization','CCC'), grid
    hold on, plot([t(1),t(end)], params.ptoTorqueOptions(:)*[1,1],'k--'), hold off, grid
    subplot(212), plot(t,thetaDot,t,thetaDotOpt_Theory), ylabel('Angular Speed [rad/s]'), xlabel('Time [s]'), legend('Numerical Optimization','|T|^2/2/R'), grid

end
end