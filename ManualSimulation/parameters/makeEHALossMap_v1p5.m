function EHA = makeEHALossMap_v1p5(A,vMax)
%%
nV = 100;
v_vals = linspace(-1,1,nV);

nT = 100;
T_vals = linspace(-1,1,nT)*3e7;

[V,T] = ndgrid(v_vals,T_vals);

Loss = NaN(size(V));
for i = 1:length(V(:))
    Q = V(i)*A;
    F = T(i)/2.7574;
    deltaP = F/A;
    Loss(i) = ehaLoss(Q,deltaP,vMax*A);
end

% Least Squares
% Reshape vectors for least squares and scale
v = reshape(V,[nV*nT,1]);
t = reshape(T,[nV*nT,1])/1e6;
l = reshape(Loss,[nV*nT,1])/1e5;

% Arrange independent variables
indep = [v.^2, v.*t, t.^2, ones(nT*nV,1)];

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
    LossHat(i) = a*V(i)^2 + b*V(i)*T(i) + c*T(i)^2 + d;
end

% output function
EHA = struct();
EHA.LossCoeffs = [a b c d];
EHA.LossFunc = @(Q,deltaP) ehaLoss(Q,deltaP,vMax*A);

if 1
    % Plot loss map
    levels = round(linspace(min(Loss(:)),max(Loss(:)),10)/1e4)*1e4/1e3;
    % figure, contour(V,T,LossHat/1e3,levels,'showtext','on'), xlabel('Speed'), ylabel('Torque'), title('Estimate')
    figure, contour(V,T,Loss/1e3,levels,'showtext','on'), xlabel('Speed'), ylabel('Torque'), title('Constant Displacement')
    % figure, surf(V,T,Loss-LossHat), xlabel('Speed'), ylabel('Torque'), title('Difference')
    % figure, surf(V,T,1e-7*T.^2), xlabel('Speed'), ylabel('Torque'), title('Electric Loss')
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

fracDisp = 1;
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*w*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
w = (Q + sign(deltaP)*QLoss) / (d*fracDisp*Scale);

T_Ideal = deltaP*d*fracDisp*Scale;
TLoss = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*w^2*rho*d^(5/3)/2)  );
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
P_L_elect = (5e-5)*T_Act^2*sign(T_Act);
P_out = w*T_Act + P_L_elect;

P_in = Q*deltaP;

Loss = abs(P_in-P_out);
if ~isfinite(Loss)
    a=1;
end
end