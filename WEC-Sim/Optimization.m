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
PR = [0 30]*1e6;
rodArea = 0.15^2*pi; % m^2: Radius squared times pi
capArea = 1.5*rodArea; % m^2: Area ratio times rod Area
r= 1.18; % Moment arm from force to torque
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
params.uDT = .1;
params.finalTime = 50;
params.t = (0:params.uDT:params.finalTime)';
params.ptoTorqueOptions = ptoTorqueOptions;

    % load ExcitationTorque.mat % Only gives forces at 0.2s intervals
    % params.T_exc = T_exc(1:length(params.t),2);
params.period = 5; % s
rampTime = 10; % s
H = freqresp(sys,2*pi/params.period);

params.sys = sys; params.r=r; params.capArea=capArea; params.rodArea=rodArea;
options = optimoptions('fmincon','PlotFcn','optimplot','MaxFunctionEvaluations',100*length(params.t),'StepTolerance',3e-5,'OptimalityTolerance',1e-7);
u0  = zeros(2,1);
A = []; b = [];
Aeq = []; beq = [];
lb = min(ptoTorqueOptions(:))*ones(size(u0))/1e6;
ub = max(ptoTorqueOptions(:))*ones(size(u0))/1e6;
nonlcon = [];
%%
if length(u0) == 1
    u_vals = linspace(min(ptoTorqueOptions(:)),max(ptoTorqueOptions(:)),25)/1e6;
    J = NaN(size(u_vals));
    Waitbar = waitbar(0,'Grid search u values');
    for uInd=1:length(u_vals)
        waitbar(uInd/length(u_vals),Waitbar)
        J(uInd) = dynamics(u_vals(uInd),params);
    end
    close(Waitbar)
    figure, plot(u_vals,J), xlabel('Control Input [Nm]'), ylabel('Cost Function [MJ]'), grid
    [~,u_opt_ind] = min(J);
    u_opt = u_vals(u_opt_ind)
elseif length(u0) == 2
    tmp = linspace(min(ptoTorqueOptions(:)),max(ptoTorqueOptions(:)),10)/1e6;
    u_vals = sort(unique([tmp,ptoTorqueOptions(:)'/1e6]));
    [u_vals1,u_vals2] = ndgrid(u_vals,u_vals); J = NaN(size(u_vals1));
    Waitbar = waitbar(0,'Grid search u values');
    for uInd1 = 1:length(u_vals)
        for uInd2 = 1:length(u_vals)
            J(uInd1,uInd2) = dynamics([u_vals1(uInd1,uInd2);u_vals2(uInd1,uInd2)],params);
        end
        waitbar(uInd1/length(u_vals),Waitbar)
    end
    close(Waitbar)
    infeasibleInds = u_vals1>0; J(infeasibleInds) = NaN;
    infeasibleInds = u_vals2<.5; J(infeasibleInds) = NaN;
    [optCost,u_opt_ind] = min(J,[],"all"); optCost
    figure, contour(u_vals1,u_vals2,J,[5,4,3,2, 1, 0, -.1, -.18 -.2]), xlabel('Control Input 1 [Nm]'), ylabel('Control Input 2 [Nm]')%, zlabel('Cost Function [MJ]'), grid
    grid, colorbar, hold on, plot([min(ptoTorqueOptions(:)),max(ptoTorqueOptions(:))]/1e6,ptoTorqueOptions(:)*[1,1]/1e6,'k--')
    plot(ptoTorqueOptions(:)*[1,1]/1e6,[min(ptoTorqueOptions(:)),max(ptoTorqueOptions(:))]/1e6,'k--')
    plot(u_vals1(u_opt_ind),u_vals2(u_opt_ind),'r*'),hold off
    % figure, surf(u_vals1,u_vals2,J), xlabel('Control Input 1 [Nm]'), ylabel('Control Input 2 [Nm]'), zlabel('Cost Function [MJ]')
    % constantplane('x',ptoTorqueOptions([1,4])/1e6,FaceAlpha=.5), constantplane('y',ptoTorqueOptions([1,4])/1e6,FaceAlpha=.5)
end
%%
params.plotFigures = 0;
u0 = zeros(50,1);

lb = min(ptoTorqueOptions(:))*ones(size(u0))/1e6;
ub = max(ptoTorqueOptions(:))*ones(size(u0))/1e6;

%lb(1) = max(lb(1),ptoTorqueOptions(1)/1e6); 
%ub(1) = min(ub(1),ptoTorqueOptions(1)/1e6);

%lb(2) = max(lb(2),ptoTorqueOptions(1)/1e6); 
%ub(2) = min(ub(1),ptoTorqueOptions(1)/1e6);

u0_vals = [min(ptoTorqueOptions(:))/1e6+.1, 0-.01, max(ptoTorqueOptions(:))/1e6-.1];
cost_vec = NaN(size(u0_vals));
uOpt_mat = NaN(length(u0),length(u0_vals));
for u0_ind = 2%1:length(u0_vals)
    u0 = u0_vals(u0_ind)*ones(size(u0));
        % uOpt = u0; costOpt = dynamics(uOpt, params);
    [uOpt,costOpt,exitflag,output] = fmincon(@(u) dynamics(u,params),u0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    cost_vec(u0_ind) = costOpt;
    uOpt_mat(:,u0_ind) = uOpt;
end
% figure, plot(uOpt_mat)
% figure, plot(u0_vals,cost_vec)
%
params.plotFigures = 1;
cost = dynamics(uOpt, params);
% uOpt

%uPI = 0*uOpt;
%cost = dynamics(uPI, params)

%
%thetaDot = eulerInt(params);
function cost = dynamics(u_scaled, params)
dt = 1e-2;
t = (0:dt:50)';
u = NaN(size(t));
u_discrete = timeseries([u_scaled;u_scaled(end)]*1e6,(0:params.uDT:params.uDT*(length(u_scaled)))+20);
warning('off','MATLAB:linearinter:noextrap'); u_zoh = resample(u_discrete,t,'zoh');


H = freqresp(params.sys,2*pi/params.period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');
rampTime = 10; % s
ramp = 0.5*(1+cos(pi+pi*t/rampTime)).*(t<rampTime) + (t>=rampTime);
Texc = 1.5e6*sin(t*2*pi/params.period).*ramp;
    % figure, plot(params.t,params.T_exc/2/Kp);
states = NaN(length(params.sys.A),length(t)); 
states(:,1) = zeros(length(params.sys.A),1);
for timeInd = 1:length(t)-1
    thetaDot(timeInd) = states(1,timeInd);
    theta(timeInd) = states(2,timeInd);
    if isfinite(u_zoh.Data(timeInd))
        u(timeInd) = u_zoh.Data(timeInd);
    else
        u(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));
    end
    states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(u(timeInd)+Texc(timeInd)))*dt;
end
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
u(timeInd+1) = -1*(Kp*thetaDot(timeInd+1) + Ki*theta(timeInd+1));
    
mechPower = u'.*thetaDot;

% compute switching losses
% Compute losses only for the segments where we're selecting the torques
capTorque = u_discrete.Data/2;% Obviously this is wrong
rodTorque = u_discrete.Data/2;
theta_discrete = resample(timeseries(theta',t),u_discrete.Time);
thetaDot_discrete = resample(timeseries(thetaDot',t),u_discrete.Time);
capLoss = SwitchingLoss(capTorque,theta_discrete.Data,thetaDot_discrete.Data,'cap',params.capArea,params.r);
rodLoss = SwitchingLoss(rodTorque,theta_discrete.Data,thetaDot_discrete.Data,'rod',params.rodArea,params.r);

cost = sum(mechPower)*dt + capLoss + rodLoss;
cost = cost/1e6;


if params.plotFigures
    Pmax = (1.5e6)^2/8/Kp
    Pact = cost/t(end)*1e6;
    Pact_kW = Pact/1e3
    thetaDotOpt_Theory = Texc/2/Kp;
    thetaOpt_Theory = cumtrapz(t,thetaDotOpt_Theory);
    figure, plot(t,u,t,-1*(Kp*thetaDotOpt_Theory + Ki*thetaOpt_Theory)), ylabel('PTO Torque [Nm]'), xlabel('Time [s]'), %legend('Numerical Optimization','CCC'), grid
    hold on, plot([t(1),t(end)], params.ptoTorqueOptions(:)*[1,1],'k--'), hold off, xlim([19 25]), grid
    % figure, plot(t,thetaDot,t,thetaDotOpt_Theory), ylabel('Angular Speed [rad/s]'), xlabel('Time [s]'), legend('Numerical Optimization','|T|^2/2/R'), grid
    % 
    % figure
    % subplot(311), plot(t,theta,t,thetaDot),grid, legend('Theta','Theta Dot'), xlabel('Time [s]'), ylabel('[rad] or [rad/s]')
    % subplot(312), plot(t,mechPower),grid, ylabel('Power [W]'), xlabel('Time [s]')
    % subplot(313), plot(t,u,t,Texc),grid, ylabel('Torque [Nm]'), xlabel('Time [s]'), legend('PTO','Excitation')
    % 
    % figure, yyaxis left,  plot(t,Texc),grid, ylabel('Exitation Torque [Nm]'), xlabel('Time [s]')
    %         yyaxis right, plot(t,thetaDot),grid, xlabel('Time [s]'), ylabel('Angular Velocity [rad/s]')
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
% a=0;b=0;c=0;d=0;

loss = sum(a.*delP.^2 + b.*Q + c.*V + d.*Q.^2)*1e3;
    % the term "a*delP^2 + b*Q + c*V + d*Q^2" is the amount of energy in a
    % 0.1 second time step because of switching.

end