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
DT = .1; % currently is both the simulation step and the optimization step
PR = [0 40]*1e6;
rodArea = 0.2^2*pi; % m^2: Radius squared times pi
capArea = 1.5*rodArea; % m^2: Area ratio times rod Area
r= 2; % distance from hinge to cylinder connection
capTorqueOptions = PR*capArea*r;
rodTorqueOptions = PR*rodArea*r;
ptoTorqueOptions = capTorqueOptions'-rodTorqueOptions;

% State space
A = [0,            -Khs/(I+Iinf),0,-1/(I+Iinf);...
    1,               0, 0,        0;...
    0,               0, 0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0, 1, -x_timeDomain(2)];
B = [1/(I+Iinf);0;0;0]; C = [1,0,0,0;0,1,0,0];
sys = ss(A,B,C,0);

params = struct;
params.plotFigures = 0;
params.DT = DT;
params.finalTime = 50;
params.t = (0:params.DT:params.finalTime)';

    % load ExcitationTorque.mat % Only gives forces at 0.2s intervals
    % params.T_exc = T_exc(1:length(params.t),2);
period = 5; % s
rampTime = 10; % s
ramp = 0.5*(1+cos(pi+pi*params.t/rampTime)).*(params.t<rampTime) + (params.t>=rampTime);
params.T_exc = 1.5e6*sin(params.t*2*pi/period).*ramp;
figure, plot(params.t,params.T_exc)
H = freqresp(sys,2*pi/period);
R = real(1/H(1));
figure, plot(params.t,params.T_exc/2/R);
Emax = (1.5e6)^2/8/R*params.finalTime

params.sys = sys; params.r=r; params.capArea=capArea; params.rodArea=rodArea;
options = optimoptions('fmincon','PlotFcn','optimplot','MaxFunctionEvaluations',200*length(params.t),'StepTolerance',1e-5,'OptimalityTolerance',1e-3);
u0  = zeros(size(params.t));
A = []; b = [];
Aeq = []; beq = [];
lb = min(ptoTorqueOptions(:))*ones(size(u0))/1e6;
ub = max(ptoTorqueOptions(:))*ones(size(u0))/1e6;
nonlcon = [];

u0_vals = [min(ptoTorqueOptions(:))/1e6+.1, 0, max(ptoTorqueOptions(:))/1e6-.1];
cost_vec = NaN(size(u0_vals));
uOpt_mat = NaN(length(params.t),length(u0_vals));
for u0_ind = 2%1:length(u0_vals)
    u0 = u0_vals(u0_ind)*ones(size(params.t));
        % uOpt = u0; costOpt = dynamics(uOpt, params);
    [uOpt,costOpt,exitflag,output] = fmincon(@(u) dynamics(u,params),u0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    cost_vec(u0_ind) = costOpt;
    uOpt_mat(:,u0_ind) = uOpt;
end
figure, plot(params.t,uOpt_mat)
figure, plot(u0_vals,cost_vec)
%%
params.plotFigures = 1;
cost = dynamics(uOpt, params)

%uPI = 0*uOpt;
%cost = dynamics(uPI, params)

%%
thetaDot = eulerInt(params);
function cost = dynamics(u_scaled, params)
u = u_scaled*1e6;

T_exc = params.T_exc;
y = lsim(params.sys,u+T_exc,params.t);
thetaDot = y(:,1);
theta = y(:,2);
    
mechPower = u.*thetaDot;

capTorque = u/2;% Obviously this is wrong
rodTorque = u/2;
capLoss = SwitchingLoss(capTorque,theta,thetaDot,'cap',params.capArea,params.r);
rodLoss = SwitchingLoss(rodTorque,theta,thetaDot,'rod',params.rodArea,params.r);

cost = sum(mechPower)*params.DT + capLoss + rodLoss;
cost = cost/1e6;

if params.plotFigures
    figure
    subplot(311), plot(params.t,theta,params.t,thetaDot),grid, legend('Theta','Theta Dot'), xlabel('Time [s]'), ylabel('[rad] or [rad/s]')
    subplot(312), plot(params.t,mechPower),grid, ylabel('Power [W]'), xlabel('Time [s]')
    subplot(313), plot(params.t,u,params.t,T_exc),grid, ylabel('Torque [Nm]'), xlabel('Time [s]'), legend('PTO','Excitation')

    figure, yyaxis left,  plot(params.t,T_exc),grid, ylabel('Exitation Torque [Nm]'), xlabel('Time [s]')
            yyaxis right, plot(params.t,thetaDot),grid, xlabel('Time [s]'), ylabel('Angular Velocity [rad/s]')

    thetaDotOpt_Theory = T_exc/2/(6.2939e7);
    thetaOpt_Theory = cumtrapz(params.t,thetaDotOpt_Theory);
    figure, plot(params.t,thetaDot,params.t,thetaDotOpt_Theory), ylabel('Angular Speed [rad/s]'), xlabel('Time [s]'), legend('Numerical Optimization','|T|^2/2/R'), grid

    H = freqresp(params.sys,2*pi/5);
    Kp = real(1/H(1)');
    Ki = -2*pi/5*imag(1/H(1)');
    figure, plot(params.t,u,params.t,-1*(Kp*thetaDotOpt_Theory + Ki*thetaOpt_Theory)), ylabel('PTO Torque [Nm]'), xlabel('Time [s]'), legend('Numerical Optimization','CCC'), grid
end
end

function loss = SwitchingLoss(u,theta,thetaDot,side,area,r)
hoseVolume = 0.1^2*pi*5;  % m^3: Radius of hose squared times pi times length of the hose
thetaMax = 45*pi/180;

% Positive thetadot results in extension of the cylinder
% Positive flow is flow into the cylinder
% Positive torque acts to extend the cylinder and move theta positive
switch side
    case 'cap'
        Q = thetaDot*r*area;
        V = (thetaMax+theta)*r*area + hoseVolume;
    case 'rod'
        Q = -thetaDot*r*area;
        V = (thetaMax-theta)*r*area + hoseVolume;
end


delP = -[u(1); diff(u)]/area/1e6;

a = 2e-4*delP.^2 - .004;
b = -0.04*delP;
c = 0;
d = 0.0060;
a=0;b=0;c=0;d=0;

loss = sum(a.*delP.^2 + b.*Q + c.*V + d.*Q.^2)*1e3;
    % the term "a*delP^2 + b*Q + c*V + d*Q^2" is the amount of energy in a
    % 0.1 second time step because of switching.

end

function thetaDot = eulerInt(params)
dt = 1e-3;
t = 0:dt:50;
period = 5; % s
rampTime = 10; % s
ramp = 0.5*(1+cos(pi+pi*t/rampTime)).*(t<rampTime) + (t>=rampTime);
Texc = 1.5e6*sin(t*2*pi/period).*ramp;
H = freqresp(params.sys,2*pi/period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');
states = NaN(length(params.sys.A),length(params.t)); u = NaN(size(t));
states(:,1) = zeros(length(params.sys.A),1);
for timeInd = 1:length(t)-1
    thetaDot(timeInd) = states(1,timeInd);
    theta(timeInd) = states(2,timeInd);
    u(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));
    states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(u(timeInd)+Texc(timeInd)))*dt;
end
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
u(timeInd+1) = -1*(Kp*thetaDot(timeInd+1) + Ki*theta(timeInd+1));
figure, plot(t,thetaDot), ylabel('Theta [rad]')
figure, plot(t,u,t,Texc), ylabel('Torque [Nm]'), legend('PTO','Excitation')

-sum(u.*thetaDot)*dt/t(end)
Pmax = (1.5e6)^2/8/Kp
tmp = 0;
end