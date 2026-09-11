% plotEHAEfficiencyMap.m
% Opt-in diagnostic plot (behind runParams.plotEHAEfficiency): contours
% the real EHA loss map (EHA.LossFunc) and compares it against the
% quadratic fit (EHA.LossCoeffs) MPC_QP actually uses internally, on the
% same (thetaDot,torque) grid convention as makeEHALossMap.m.
% Calls: none
% Called by: parameters/getHydraulic.m
function plotEHAEfficiencyMap(EHA, capArea)
% Visualizes the EHA loss map from makeEHALossMap.m: the true pump-physics
% loss (via EHA.LossFunc, the fixed-2000RPM variable-displacement model)
% over a (speed,torque) grid, and a second figure comparing it against the
% quadratic fit (EHA.LossCoeffs) that the electrical-optimized MPC_QP
% controller actually uses internally (ctrl.MPC in getControl.m).
%
% The (speed W [rad/s], torque T [Nm]) -> (flow Q, deltaP) mapping here
% matches makeEHALossMap.m's own internal grid exactly (V=W*2.7574,
% Q=V*capArea, F=T/2.7574, deltaP=F/capArea) so this lines up with what the
% controller/evaluation code actually sees.

nW = 60; W_vals = linspace(-.5,.5,nW);
nT = 60; T_vals = linspace(-1,1,nT)*3e7;
[W,T] = ndgrid(W_vals,T_vals);

a = EHA.LossCoeffs(1); b = EHA.LossCoeffs(2); c = EHA.LossCoeffs(3); d = EHA.LossCoeffs(4);

Loss = NaN(size(W));
LossHat = NaN(size(W));
for i = 1:numel(W)
    V = W(i)*2.7574;
    Q = V*capArea;
    F = T(i)/2.7574;
    deltaP = F/capArea;
    Loss(i) = EHA.LossFunc(Q,deltaP);
    LossHat(i) = a*W(i)^2 + b*W(i)*T(i) + c*T(i)^2 + d;
end

levels = round(linspace(min(Loss(:)),max([Loss(:);LossHat(:)]),12)/1e3)*1e3/1e3;

% Figure 1: the efficiency/loss map on its own
figure('Name','EHA Loss Map');
contourf(W,T/1e6,Loss/1e3,levels,'ShowText','on');
colorbar
xlabel('Equivalent shaft speed, \theta_{dot} [rad/s]')
ylabel('Torque [MNm]')
title({'EHA Loss Map','(pump physics, 2000 RPM constant speed)','[kW]'})

% Figure 2: map vs. quadratic-fit comparison
figure('Name','EHA Loss Map vs Quadratic Fit');
subplot(1,2,1)
contourf(W,T/1e6,Loss/1e3,levels,'ShowText','on');
colorbar, xlabel('\theta_{dot} [rad/s]'), ylabel('Torque [MNm]')
title('Real (pump physics)')
subplot(1,2,2)
contourf(W,T/1e6,LossHat/1e3,levels,'ShowText','on');
colorbar, xlabel('\theta_{dot} [rad/s]'), ylabel('Torque [MNm]')
title('Quadratic fit (used by MPC\_QP controller)')
sgtitle('EHA Loss: Real vs Quadratic Fit [kW]')

end
