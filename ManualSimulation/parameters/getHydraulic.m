function hyd = getHydraulic(~)

% Set Valued Control Inputs
hyd.pressureRails = [0 30]*1e6;
hyd.stroke = 5;
hyd.rodArea = 0.2^2*pi; % m^2: Radius squared times pi
hyd.capArea = 1.5*hyd.rodArea; % m^2: Area ratio times rod Area

% Congifuration of hydraulic cylinder
hyd.r= 1.18; % Moment arm from force to torque

cylHorizDist = 7; % Distance from hinge to cylinder along sea floor
r_cyl = 3; % distance from hinge to cylinder along flap
hyd.L = @(theta) sqrt(cylHorizDist^2+r_cyl^2 + 2*cylHorizDist*r_cyl*sin(theta));
hyd.dLdt = @(theta,thetaDot) cylHorizDist*r_cyl*cos(theta).*thetaDot./hyd.L(theta);

hyd.L_equilib = hyd.L(0);
hyd.L_retract = hyd.L_equilib - hyd.stroke/2; % Length of cylinder at full retraction

capTorqueOptions = hyd.pressureRails * hyd.capArea * hyd.r;
rodTorqueOptions = hyd.pressureRails * hyd.rodArea * hyd.r;
hyd.ptoTorqueOptions = capTorqueOptions'-rodTorqueOptions;

end