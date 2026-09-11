% makeEHALossMap.m
% Baseline EHA model: constant speed (2000 RPM), variable displacement.
% Local ehaLoss solves for fracDisp to match commanded flow, computes
% torque loss (viscous/Coulomb/windage) and i^2R copper loss. Fits the
% quadratic loss surrogate (LossCoeffs) against (thetaDot [rad/s], torque
% [Nm]) directly, matching how MPC_EHA (getControl.m) uses it.
% Calls: none
% Called by: parameters/getHydraulic.m, plotEHAEfficiencyMap.m (for the
%   real-vs-fit comparison plot), makeEHALossMap_fixedDisp.m (references
%   its ehaLoss physics in a comment only, does not call it)
function EHA = makeEHALossMap(A,vMax)
%%
nw = 100;
w_vals = linspace(-.5,.5,nw);

nT = 100;
T_vals = linspace(-1,1,nT)*3e7;

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
% Reshape vectors for least squares and scale
w = reshape(W,[nw*nT,1]);
t = reshape(T,[nw*nT,1])/1e6;
l = reshape(Loss,[nw*nT,1])/1e5;

% Arrange independent variables
indep = [w.^2, w.*t, t.^2, ones(nT*nw,1)];

% Solve for the lease squares
par = pinv(indep)*l;

% rescale the coefficients
coeffs = par*1e5;
a = coeffs(1);
b = coeffs(2)/1e6;
c = coeffs(3)/1e6/1e6;
d = coeffs(4);

LossHat = NaN(size(Loss));
for i = 1:length(Loss(:))
    LossHat(i) = a*W(i)^2 + b*W(i)*T(i) + c*T(i)^2 + d;
end

% output function
EHA = struct();
EHA.LossCoeffs = [a b c d];
EHA.LossFunc = @(Q,deltaP) ehaLoss(Q,deltaP,vMax*A);

if 0
    % Plot loss map
    levels = round(linspace(min(Loss(:)),max(Loss(:)),10)/1e4)*1e4/1e3;
    figure, contour(W,T,LossHat/1e3,levels,'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Estimate')
    figure, contour(W,T,Loss/1e3,levels,'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Constant Speed')
    % figure, surf(W,T,Loss-LossHat), xlabel('Speed'), ylabel('Torque'), title('Difference')
    % figure, surf(W,T,1e-7*T.^2), xlabel('Speed'), ylabel('Torque'), title('Electric Loss')
end
a=1;
end

function Loss = ehaLoss(Q,deltaP,maxFlow)
% Define Pump Constants
% Angular Velocity
Wrpm = 2000; %revolutions per minute
w = Wrpm.*(2*pi/60); % radians per second

Scale = maxFlow/w*2*pi*1e6/107 *1.2; % The 1.2 is to account for losses (oversize so we can actually hit the max flow at 2000 RPM)

% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
    D = 107; % cc/rev
    d = (D*100^-3)/(2*pi); % m^3/rad 
% Torques Loss Constants
    Cf =  53.7e-3;
    Ch = 53.6;
    Cv = 23.5e3;
% Flow Loss Constants
    Cs = 4.26e-9;
    Cst = 0*1e-5;

mu=(32e-6)*870;
B = 1.7e9;
rho = 870;

% Calculate fracDisp assuming it is a posative value
fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale-sign(deltaP)*Scale*abs(d*w*deltaP/B));
if fracDisp <= 0 % if the assumption that fracdisp is + is incorrect, recalculate fracDisp assuming fracDisp is -
    fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale+sign(deltaP)*Scale*abs(d*w*deltaP/B));
end

% Check if fracDisp was calculated correctly
% Calculate the error due to the fractional displacement calculation
% If I use the fracDisp I calculated to solve for Q_Act (which was given)
% Do I get the same values?
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*w*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
Q_Act_calc = w*d*fracDisp*Scale - sign(deltaP)*QLoss;

T_Ideal = deltaP*d*fracDisp*Scale;
TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );

Q_Ideal = fracDisp*d*w*Scale;
T_Act = T_Ideal + sign(w)*TLoss; % |T_Act| needs be < |T_Ideal|

% Power out with i^r losses (copper loss, proportional to T_Act^2)
P_L_elect = (5e-5)*T_Act^2*sign(T_Act);
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;

Loss = abs(P_in-P_out);
if ~isfinite(Loss)
    a=1;
end
end