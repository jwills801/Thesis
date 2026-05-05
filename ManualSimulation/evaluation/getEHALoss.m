function eval = getEHALoss(params,dyn)
%% Calculate volume and flow in each side
[cap,~] = params.hyd.getVolandFlow(params,dyn.states);
Q = cap.velA;

theta = dyn.states(2,:)';
thetaDot = dyn.states(1,:)';
F = dyn.u ./ params.hyd.Force2Torque(theta);
deltaP = F/params.hyd.capArea;
vel = params.hyd.dLdt(theta,thetaDot);

loss = NaN(length(dyn.t),1);
for t_ind = 1:length(dyn.t)
    loss(t_ind) = params.hyd.EHA.LossFunc(Q(t_ind),deltaP(t_ind));
    
end
lossAfterRamp = loss(dyn.t > params.simu.rampTime);

% Calculate loss from the convex approximation
a = params.hyd.EHA.LossCoeffs(1); b = params.hyd.EHA.LossCoeffs(2);
c = params.hyd.EHA.LossCoeffs(3); d = params.hyd.EHA.LossCoeffs(4);
lossHat = a*vel.^2 + b*vel.*dyn.u + c*dyn.u.^2 + d;
lossHatAfterRamp = lossHat(dyn.t > params.simu.rampTime);
% figure, plot(dyn.t,loss,dyn.t,lossHat), legend('Real','Approximation')

% Output results
eval.TotalLoss = sum(loss)*params.simu.dt;
eval.loss = loss;
eval.TotalLossAfterRamp = sum(lossAfterRamp)*params.simu.dt;
eval.aveLoss = eval.TotalLossAfterRamp / (params.simu.finalTime - params.simu.rampTime);
eval.TotalLossHatAfterRamp = sum(lossHatAfterRamp)*params.simu.dt;
eval.aveLossHat = eval.TotalLossHatAfterRamp / (params.simu.finalTime - params.simu.rampTime);
end