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
params.hoseVolume = params.capArea*5; % [m^3] 5 meters of hose the same diameter as the cap


%% For active only
% Kp = 7.7426e+06; Ki = 7.7426e+06; % T8 H2.5
Kp = 2.1544e6; Ki = 1;
valveSettlingTime = .1; % [s]
valveTransferFunction = tf(4/valveSettlingTime,[1 4/valveSettlingTime]);


%% For passive only
params.pressure = 7e6;

%% for HHEA only
load AP_107cc_Sept16.mat
params.Disp = 3000*1e-6/2/pi; 
params.scaleHECM = params.Disp/Disp; % The map assumes a 107cc motor
params.P1_Mapping = P1_Mapping;
params.w1rad_Mapping = w1rad_Mapping;
params.T1_Act_Mapping = T1_Act_Mapping;
params.Q1_Act_Mapping = Q1_Act_Mapping;
params.shaftInertia = 1/2*10*.01^2; % [kg m^2] 10 kg, 0.05m radius motor shaft