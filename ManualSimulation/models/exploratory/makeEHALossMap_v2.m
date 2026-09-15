% makeEHALossMap_v2.m
% Exploratory EHA model: both speed and displacement free -- local
% ehaLoss sweeps 1000 candidate shaft speeds per (Q,deltaP) point and
% picks whichever (speed,displacement) combination minimizes loss
% (subject to |fracDisp|<=1). An idealized best-case bound, not a
% controller-ready model (no LossCoeffs unit issue since it's never fed
% to MPC_EHA). Uses a 1.5x pump-oversizing margin vs. the 1.2x in the
% other two files -- a known inconsistency, not reconciled.
% Calls: none
% Called by: none currently wired into a live path (used ad hoc from
%   scratch comparison scripts during development)
function EHA = makeEHALossMap_v2(A,vMax)
%%
nthetaDot = 50;
thetaDot_vals = linspace(-.4,.4,nthetaDot);

nT = 100;
T_vals = linspace(-3,3,nT)*1e7;

[ThetaDot,T] = ndgrid(thetaDot_vals,T_vals);

Loss = NaN(size(ThetaDot));
Chi = NaN(size(ThetaDot));
Speed = NaN(size(ThetaDot));
for i = 1:length(ThetaDot(:))
    V = ThetaDot(i)*2.7574;
    Q = V*A;
    F = T(i)/2.7574;
    deltaP = F/A;
    [Loss(i),Chi(i),Speed(i)] = ehaLoss(Q,deltaP,vMax*A);
end

% Round Chi if close to one
Chi(Chi>.99)=1;
Chi(Chi<-.99)=-1;

% Least Squares
% Reshape vectors for least squares and scale
v = reshape(ThetaDot,[nthetaDot*nT,1]);
t = reshape(T,[nthetaDot*nT,1])/1e6;
l = reshape(Loss,[nthetaDot*nT,1])/1e5;

% Arrange independent variables
indep = [v.^2, v.*t, t.^2, ones(nT*nthetaDot,1)];

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
    LossHat(i) = a*ThetaDot(i)^2 + b*ThetaDot(i)*T(i) + c*T(i)^2 + d;
end

% output function
EHA = struct();
EHA.LossCoeffs = [a b c d];
EHA.LossFunc = @(Q,deltaP) ehaLoss(Q,deltaP,vMax*A);

if 0
    % Plot loss map
    levels = round(linspace(min(Loss(:)),max(Loss(:)),10)/1e4)*1e4/1e3;
    figure, contour(ThetaDot,T,LossHat/1e3,levels,'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Estimate')
    figure, contour(ThetaDot,T,Loss/1e3,levels,'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Real Deal')

    figure, contour(ThetaDot,T,Chi,[-.95 -.8 -.1 .1 .8 .95],'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Fractional Disp')
    figure, surf(ThetaDot,T,Chi), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Fractional Disp')
    figure, contour(ThetaDot,T,Speed*60/2/pi,'showtext','on'), xlabel('Speed [rad/s]'), ylabel('Torque'), title('Shaft Speed [RPM]')
    % figure, surf(V,T,Loss-LossHat), xlabel('Speed'), ylabel('Torque'), title('Difference')
    % figure, surf(V,T,1e-7*T.^2), xlabel('Speed'), ylabel('Torque'), title('Electric Loss')
end

end

function [Loss,Chi,Speed] = ehaLoss(Q,deltaP,maxFlow)
% Define Pump Constants
% Angular Velocity
Wmaxrpm = 2000; %revolutions per minute
wmax = Wmaxrpm.*(2*pi/60); % radians per second

Scale = maxFlow/wmax*2*pi*1e6/107 *1.5; % The 1.2 is to account for losses (oversize so we can actually hit the max flow at 2000 RPM)

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

wVals = linspace(10,wmax,1000);
LossVec = NaN(size(wVals));
for wInd = 1:length(wVals)
    w = wVals(wInd);
    % Calculate fracDisp assuming it is a posative value
    fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale-sign(deltaP)*Scale*abs(d*w*deltaP/B));
    if fracDisp <= 0 % if the assumption that fracdisp is + is incorrect, recalculate fracDisp assuming fracDisp is -
        fracDisp = (Q + sign(deltaP)*Scale*abs(d*Cs*(deltaP)/mu) + sign(deltaP)*Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5))/(w*d*Scale+sign(deltaP)*Scale*abs(d*w*deltaP/B));
    end
    fracDispVec(wInd) = fracDisp;

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

    % Power out with i^r losses
    P_L_elect = (5e-5)*T_Act^2*sign(T_Act);
    P_out = w*T_Act + P_L_elect;

    P_in = Q*deltaP;
    LossVec(wInd) = abs(P_in-P_out);
    if abs(fracDisp) > 1
        LossVec(wInd) = NaN;
    end
end
[Loss,I] = min(LossVec);
Chi = fracDispVec(I);
Speed = wVals(I);

if ~isfinite(Loss)
    a=1;
end
end