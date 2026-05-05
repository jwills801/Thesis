A = .2;

nV = 100;
v_vals = linspace(-1,1,nV);

nT = 100;
T_vals = linspace(-1,1,nT)*1e6;

[V,T] = ndgrid(v_vals,T_vals);

Loss = NaN(size(V));
for i = 1:length(V(:))
    Q = V(i)*A;
    F = T(i)/2.7574;
    deltaP = F/A;
    Loss(i) = ehaLoss(Q,deltaP,.5*A);
end

% Least Squares
v = reshape(V,[nV*nT,1]);
t = reshape(T,[nV*nT,1])/1e6;
l = reshape(Loss,[nV*nT,1])/1e5;

indep = [v.^2, 0*v.*t, t.^2, ones(nT*nV,1)];
par = pinv(indep)*l;
a = par(1);
b = par(2);
c = par(3);
d = par(4);

LossHat = NaN(size(Loss));
for i = 1:length(Loss(:))
    T_scaled = T(i)/1e6;
    L_unscaled = a*V(i)^2 + 0*b*V(i)*T_scaled + c*T_scaled^2 + d;
    LossHat(i) = L_unscaled*1e5;
end

% output function
EHA.LossCoeffs = par;
EHA.LossFunc = @(Q,deltaP) ehaLoss(Q,deltaP,1*A);

if 0
    % Plot loss map
    levels = round(linspace(min(Loss(:)),max(Loss(:)),10)/1e4)*1e4/1e3;
    figure, contour(V,T,LossHat/1e3,levels,'showtext','on'), xlabel('Speed'), ylabel('Torque'), title('Estimate')
    figure, contour(V,T,Loss/1e3,levels,'showtext','on'), xlabel('Speed'), ylabel('Torque'), title('Real Deal')
    figure, surf(V,T,Loss-LossHat), xlabel('Speed'), ylabel('Torque'), title('Difference')
    figure, surf(V,T,1e-7*T.^2), xlabel('Speed'), ylabel('Torque'), title('Electric Loss')
end



function Loss = ehaLoss(Q,deltaP,maxFlow)
% Define Pump Constants
% Angular Velocity
Wrpm = 2000; %revolutions per minute
w = Wrpm.*(2*pi/60); % radians per second

Scale =maxFlow/w*2*pi*1e6/107;

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

% Power out with 90% effiency
if T_Act < 0
    P_out = .9*w*T_Act;
    P_L_elect = -.1*w*T_Act; % the negative sign makes sure the loss is posative
else
    P_out = w*T_Act/.9;
    P_L_elect =  (1/.9-1)*w*T_Act; % This loss accounts for energy that needs to come FROM the generator to the system.
end

% Power out with i^r losses
P_L_elect = (1e-7)*T_Act^2*sign(T_Act);
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;

Loss = abs(P_in-P_out);
if ~isfinite(Loss)
    a=1;
end
end