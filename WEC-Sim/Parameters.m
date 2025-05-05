%% For all PTOS
params = struct();
params.capArea = (50*.025)^2/4;
     % This is a four foot diameter bore. That's huge.
     % This could be split up into multiple actuators,
     % And/or we could lower the connection point to the flap
params.rodArea = params.capArea/1.5;
params.cylinderStroke = cylinderStroke;
params.beta = 1.8e9;
params.valveConstant = 1e-3;

%% For active only
Kp = 7.7426e+06; Ki = 7.7426e+06;

%% For passive only
params.pressure = 7e6;