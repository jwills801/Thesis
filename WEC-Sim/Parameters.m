params = struct();
params.capArea = (50*.025)^2/4;
     % This is a four foot diameter bore. That's huge.
     % This could be split up into multiple actuators,
     % And/or we could lower the connection point to the flap
params.rodArea = params.capArea/1.5;
params.cylinderStroke = 5;
params.beta = 1.8e9;
params.valveConstant = 1e-4;

%params.pressure = pressure;
% Kp = 1:0.5:5; Ki = 5:1:15;