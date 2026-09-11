% makeEHALossMap_fixedDisp.m
% Controller-ready fixed-displacement (chi=1), variable-speed EHA model
% -- same physics as makeEHALossMap_v1p5.m's ehaLoss, but LossCoeffs are
% fit against (thetaDot,torque) directly (matching makeEHALossMap.m's
% convention), avoiding v1p5's unit mismatch. Used together with
% getPhysical.m's reflected-inertia term (both keyed off
% runParams.ehaFixedDisplacement).
% Calls: none
% Called by: parameters/getHydraulic.m
function EHA = makeEHALossMap_fixedDisp(A,vMax)
% Fixed-displacement (chi=1), variable-speed EHA loss map. Shaft speed is
% solved (not fixed at 2000 RPM) to match the commanded flow at fracDisp=1
% -- same physics as makeEHALossMap_v1p5.m's ehaLoss, but LossCoeffs here
% are fit against (thetaDot [rad/s], torque [Nm]) directly, matching what
% MPC_EHA (control/getControl.m) actually multiplies them against. (v1p5
% fits against cylinder velocity V=2.7574*thetaDot instead -- a unit
% mismatch if plugged into the controller as-is, off by a factor of
% 2.7574^2 on the quadratic-in-speed term. This file avoids that.)
nw = 100; w_vals = linspace(-.5,.5,nw);

nT = 100; T_vals = linspace(-1,1,nT)*3e7;

[W,T] = ndgrid(w_vals,T_vals);

Loss = NaN(size(W));
for i = 1:length(W(:))
    V = W(i)*2.7574;
    Q = V*A;
    F = T(i)/2.7574;
    deltaP = F/A;
    Loss(i) = ehaLoss(Q,deltaP,vMax*A);
end

% Least Squares
w = reshape(W,[nw*nT,1]);
t = reshape(T,[nw*nT,1])/1e6;
l = reshape(Loss,[nw*nT,1])/1e5;

indep = [w.^2, w.*t, t.^2, ones(nT*nw,1)];
par = pinv(indep)*l;

coeffs = par*1e5;
a = coeffs(1);
b = coeffs(2)/1e6;
c = coeffs(3)/1e6/1e6;
d = coeffs(4);

EHA = struct();
EHA.LossCoeffs = [a b c d];
EHA.LossFunc = @(Q,deltaP) ehaLoss(Q,deltaP,vMax*A);
end

function Loss = ehaLoss(Q,deltaP,maxFlow)
% Same pump physics/coefficients as makeEHALossMap.m and _v1p5.m, but
% fracDisp fixed at 1 (max displacement) and shaft speed w solved to match
% the commanded flow -- i.e. the quasi-static "fixed displacement,
% variable speed" pump, evaluated at whatever (Q,deltaP) the fitting grid
% or a real sim trajectory hands it.
Wrpm = 2000; % nominal/rated speed used ONLY to size the pump (Scale), not the actual operating speed
w0 = Wrpm.*(2*pi/60);

Scale = maxFlow/w0*2*pi*1e6/107 *1.2;

D = 107; % cc/rev
d = (D*100^-3)/(2*pi);
Cf =  53.7e-3;
Ch = 53.6;
Cv = 23.5e3;
Cs = 4.26e-9;
Cst = 0*1e-5;

mu=(32e-6)*870;
B = 1.7e9;
rho = 870;

fracDisp = 1;
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*w0*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
w = (Q + sign(deltaP)*QLoss) / (d*fracDisp*Scale);

TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );
T_Ideal = deltaP*d*fracDisp*Scale;
T_Act = T_Ideal + sign(w)*TLoss;

% Copper loss (i^2*R, proportional to T_Act^2) -- see makeEHALossMap.m
P_L_elect = (5e-5)*T_Act^2*sign(T_Act);
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;
Loss = abs(P_in-P_out);
if ~isfinite(Loss)
    a=1;
end
end
