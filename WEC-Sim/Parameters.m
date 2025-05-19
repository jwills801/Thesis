function params = Parameters(~)
%% For all PTOS
params = struct();
params.capArea = (50*.025)^2/4;
     % This is a four foot diameter bore. That's huge.
     % This could be split up into multiple actuators,
     % And/or we could lower the connection point to the flap
params.rodArea = params.capArea/1.5;
params.cylinderStroke = 3;
params.beta = 1.8e9;
params.hoseVolume = params.capArea*5; % [m^3] 5 meters of hose the same diameter as the cap - possibly not used

maxActuatorSpeed = 1; % [m/s]
maxFlow = maxActuatorSpeed*params.capArea;
acceptablePressureDrop = 1e6;
params.valveConstant = maxFlow/sqrt(acceptablePressureDrop);


%% For active only
% Control gains
    % Simple PTO blocks just use Kp and Ki - not in the params structure - because it makes the mcr easier.
    % Kp = 7.7426e+06; Ki = 7.7426e+06; % T8 H2.5
    % params.Kp = 2.1544e6; params.Ki = 1; % T20 H2.5
    params.Kp = 4e6; params.Ki = -1e6; % T20 H4.2

valveSettlingTime = .1; % [s]
valveTransferFunction = tf(4/valveSettlingTime,[1 4/valveSettlingTime]);
params.valveNumerator = valveTransferFunction.Numerator{1};
params.valveDenominator = valveTransferFunction.Denominator{1};


%% For passive only
params.pressure = 7.5e6; % From grid search close to resonance
params.checkValveCrackingPressure = 1e5;
params.checkValveStroke = .0254; % 1 in stroke


%% for HHEA only
load AP_107cc_Sept16.mat
params.maxMotorSpeed = 200; % [rad/s]
% calculate required HECM size
params.Disp = maxActuatorSpeed*params.rodArea/params.maxMotorSpeed; % [m^3/rad]

params.scaleHECM = params.Disp/Disp; % The map assumes a 107cc motor
params.P1_Mapping = P1_Mapping; params.w1rad_Mapping = w1rad_Mapping; params.T1_Act_Mapping = T1_Act_Mapping; params.Q1_Act_Mapping = Q1_Act_Mapping;
params.P_Map = P_Map; params.Q_Map = Q_Map; params.W_Map = W_Map;
%params.shaftInertia = 1/2*10*.01^2; % [kg m^2] 10 kg, 0.05m radius motor shaft
end