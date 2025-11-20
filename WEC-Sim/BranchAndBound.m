clear, close all
params = getParameters;
U = params.PR'/1e6;
unsolvedNodes = {};
solvedNodes = {};
timeSteps = 6;
unsolvedNodes{1}.lb = min(U)*ones(timeSteps*2,1);
unsolvedNodes{1}.ub = max(U)*ones(timeSteps*2,1);
unsolvedNodes{1}.parentCost = NaN;
bestCost = 1e3;
flag = 0;
iter = 1;
iterMax = 1000;
%%
tic
while flag == 0
% take next node to solve
    node = unsolvedNodes{1};
toc
% solve node
tic
    node = solveNode(node,params);
toc

% Remove solved node
    unsolvedNodes(1) = [];
    solvedNodes{end+1} = node;

% At the beginning, round to the nearest node to get an integer solution
if iter == 1
    [M,I] = min(abs(solvedNodes{1}.u' - repmat(U,[1,timeSteps*2])));
    quickNode.lb = U(I);
    quickNode.ub = U(I);

    quickNode.lb = max(U)*ones(timeSteps*2,1);
    quickNode.ub = max(U)*ones(timeSteps*2,1);

    quickNode.parentCost = solvedNodes{1}.cost;
    quickNode = solveNode(quickNode,params);
    bestCost = quickNode.cost;
    solvedNodes{end+1} = quickNode;
    bestNodeInd = length(solvedNodes);
end

% Prune or branch
    if node.cost > bestCost % prune
        disp(['Pruned at node ', num2str(length(solvedNodes))])
    else % branch
        % Find out which variables are already in the allowable set
        node.inAllowableSet = min(abs(node.u' - repmat(U,[1,timeSteps*2]))) < .001;

        % Branch on the first variable not in the allowable set
        branchingVariable = find(~node.inAllowableSet,1);
        if isempty(branchingVariable) % Solution that is completely in the allowable set
            % This is the best cost. If it wasn't we wouldn't be in the
            % branch part of this if statement
            bestCost = node.cost;
            bestNodeInd = length(solvedNodes);

            % Prune nodes whose parent are already worse than this new best
                % Children nodes will always be worse than their parents
            if ~isempty(unsolvedNodes)
                tmp = [unsolvedNodes{:}];
                unsolvedNodes([tmp.parentCost] > bestCost) = [];
            end
        else
        newNodeHigh.ub = node.ub;
        newNodeHigh.lb = node.lb;
        newNodeLow = newNodeHigh;
    % find which index to branch
        tmp = node.u(branchingVariable) - U;
        I = find(tmp>0,1,'last');
        % I and (I+1) are the indices of the discrete input set that we're
        % bounding on
        newNodeHigh.lb(branchingVariable) = U(I+1);
        newNodeLow.ub(branchingVariable) = U(I);

        newNodeHigh.parentCost = node.cost;
        newNodeLow.parentCost = node.cost;

        unsolvedNodes{end+1} = newNodeHigh;
        unsolvedNodes{end+1} = newNodeLow;

        end
    end

iter = iter + 1;
disp(['Nodes Complete: ',num2str(length(solvedNodes))])
disp(['Current Unsolved Nodes: ',num2str(length(unsolvedNodes))])

% Check to see if we should stop
    if isempty(unsolvedNodes)
        flag = 1;
        disp('Solution found')
        bestNode = solvedNodes{bestNodeInd}
    elseif iter > iterMax
        flag = -1;
        disp('Max iterations')
    end
end
toc
params.plotFigures = 1;
%%
tic
cost = dynamics(solvedNodes{bestNodeInd}.u,params)
toc
return
%%
params = getParameters;
U = params.PR'/1e6;
node = struct;
timeSteps = 2;
node.lb = min(U)*ones(timeSteps*2,1);
node.ub = max(U)*ones(timeSteps*2,1);

node = solveNode(node,params)

tic
params.plotFigures = 1;
cost = dynamics(node.u,params);
toc

function node = solveNode(node,params)
% options = optimoptions('fmincon','PlotFcn','optimplot',OptimalityTolerance=1e-8);
options = optimoptions('fmincon');
A = []; b = [];
Aeq = []; beq = [];
lb = node.lb;
ub = node.ub;
nonlcon = [];
u0=0*ones(size(lb));

[uOpt,costOpt,exitflag,output] = fmincon(@(u) dynamics(u,params),u0,A,b,Aeq,beq,lb,ub,nonlcon,options);
node.u = uOpt;
node.cost = costOpt;

costOpt
variables = [lb uOpt ub]

        params.plotFigures = 1;
        cost = dynamics(node.u,params);
        drawnow
        params.plotFigures = 0;
end

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

function cost = dynamics(u_scaled,params)
dt = 1e-2;
controlStartTime = 25;
t = (0:dt:params.finalTime)';
TorquePTO = NaN(size(t));
nU = floor(length(u_scaled)/2);

% u_scaled_cap = u_scaled(1:nU);
% u_scaled_rod = u_scaled((nU+1:end));
u_scaled_cap = u_scaled(1:2:end-1);
u_scaled_rod = u_scaled(2:2:end);

u_discrete_cap = timeseries([u_scaled_cap;u_scaled_cap(end)]*1e6,(0:params.uDT:(params.uDT*nU))+controlStartTime);
u_discrete_rod = timeseries([u_scaled_rod;u_scaled_rod(end)]*1e6,(0:params.uDT:(params.uDT*nU))+controlStartTime);
warning('off','MATLAB:linearinter:noextrap'); uCap_zoh = resample(u_discrete_cap,t,'zoh'); uRod_zoh = resample(u_discrete_rod,t,'zoh');

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
    if isfinite(uCap_zoh.Data(timeInd))
        TorquePTO(timeInd) = uCap_zoh.Data(timeInd)*params.capArea*params.r - uRod_zoh.Data(timeInd)*params.rodArea*params.r;
    else
        TorquePTO(timeInd) = -1*(Kp*thetaDot(timeInd) + Ki*theta(timeInd));
    end
    states(:,timeInd + 1) = states(:,timeInd) + (params.sys.A*states(:,timeInd) + params.sys.B*(TorquePTO(timeInd)+Texc(timeInd)))*dt;
end
thetaDot(timeInd+1) = states(1,timeInd+1);
theta(timeInd+1) = states(2,timeInd+1);
TorquePTO(timeInd+1) = -1*(Kp*thetaDot(timeInd+1) + Ki*theta(timeInd+1));
    
mechPower = TorquePTO'.*thetaDot;

% compute switching losses
% Compute losses only for the segments where we're selecting the torques
theta_discrete = resample(timeseries(theta',t),u_discrete_cap.Time);
thetaDot_discrete = resample(timeseries(thetaDot',t),u_discrete_cap.Time);
capLoss = SwitchingLoss(u_discrete_cap.Data,theta_discrete.Data,thetaDot_discrete.Data,'cap',params.capArea,params.r);
rodLoss = SwitchingLoss(u_discrete_rod.Data,theta_discrete.Data,thetaDot_discrete.Data,'rod',params.rodArea,params.r);

cost = sum(mechPower)*dt + capLoss + rodLoss;

cost = cost/1e6;


if params.plotFigures
    Pmax = (1.5e6)^2/8/Kp;
    Pact = cost/t(end)*1e6;
    Pact_kW = Pact/1e3;
    thetaDotOpt_Theory = Texc/2/Kp;
    thetaOpt_Theory = cumtrapz(t,thetaDotOpt_Theory);
    figure, subplot(211), plot(t,TorquePTO,t,-1*(Kp*thetaDotOpt_Theory + Ki*thetaOpt_Theory)), ylabel('PTO Torque [Nm]'), xlabel('Time [s]'), %legend('Numerical Optimization','CCC'), grid
    hold on, plot([t(1),t(end)], params.ptoTorqueOptions(:)*[1,1],'k--'), hold off, xlim([controlStartTime-1 controlStartTime+5]), grid
    subplot(212), plot(t,thetaDot,t,thetaDotOpt_Theory), ylabel('Angular Speed [rad/s]'), xlabel('Time [s]'), legend('Numerical Optimization','|T|^2/2/R'), xlim([controlStartTime-1 controlStartTime+5]), grid

    % figure
    % subplot(311), plot(t,theta,t,thetaDot),grid, legend('Theta','Theta Dot'), xlabel('Time [s]'), ylabel('[rad] or [rad/s]')
    % subplot(312), plot(t,mechPower),grid, ylabel('Power [W]'), xlabel('Time [s]')
    % subplot(313), plot(t,TorquePTO,t,Texc),grid, ylabel('Torque [Nm]'), xlabel('Time [s]'), legend('PTO','Excitation')
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

delP = -[0;diff(u)]/area/1e6;

a = 2e-4*delP.^2 - .004;
b = -0.04*delP;
c = 0;
d = 0.0060;
% a=0;b=0;c=0;d=0;

loss = sum(a.*delP.^2 + b.*Q + c.*V + d.*Q.^2);
    % the term "a*delP^2 + b*Q + c*V + d*Q^2" is the amount of energy in a
    % 0.1 second time step because of switching.

end