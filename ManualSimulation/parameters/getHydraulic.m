function hyd = getHydraulic(~)

% Set Valued Control Inputs
hyd.pressureRails = [0 30]*1e6;
hyd.stroke = 5;
hyd.rodArea = 0.2^2*pi; % m^2: Radius squared times pi
hyd.capArea = 1.5*hyd.rodArea; % m^2: Area ratio times rod Area

hyd.r= 1.18; % Moment arm from force to torque

capTorqueOptions = hyd.pressureRails * hyd.capArea * hyd.r;
rodTorqueOptions = hyd.pressureRails * hyd.rodArea * hyd.r;
hyd.ptoTorqueOptions = capTorqueOptions'-rodTorqueOptions;

end