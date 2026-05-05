function eval = getEHALoss(params,dyn)
% Calculate volume and flow in each side
[cap,~] = params.hyd.getVolandFlow(params,dyn.states);
Q = cap.velA;

theta = dyn.states(2,:);
F = dyn.u ./ params.hyd.Force2Torque(theta);
deltaP = F/params.hyd.capArea;

loss = NaN(length(dyn.t),1);
for t_ind = 1:length(dyn.t)
    loss(t_ind) = params.hyd.EHA.LossFunc(Q(t_ind),deltaP(t_ind));
end
lossAfterRamp = loss(dyn.t > params.simu.rampTime);

% Output results
eval.TotalEHALoss = sum(loss)*params.simu.dt;
eval.loss = loss;
eval.TotalEHALossAfterRamp = sum(lossAfterRamp)*params.simu.dt;
eval.aveEHALoss = eval.TotalEHALossAfterRamp / (params.simu.finalTime - params.simu.rampTime);
% params.hyd.EHA.LossCoeffs
end