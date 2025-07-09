% Currently configured for debugging F=Ma
    % No PTO force, no waves, just an initial condition.

peakWavePeriod = 20;
codesign = 0;
% Simple force laws
    % PTO = 'Continuous PI'; Kp = 4e6; Ki = 3e6;
    % PTO = 'Discrete PI'; Kp = 4e6; Ki = 3e6;
    % PTO = 'Rectifying'; % pressure = 5e6;

    % Turn off PTO force for debugging
    PTO = 'Continuous PI'; Kp = 0; Ki = 0;

    Kp = 5e5; Ki = -1e5;

% Hydraulic PTO models
    % PTO = 'Active Valving';
    % PTO = 'EHA'; codesign = 1;
    % PTO = 'HHEA';
    % PTO = 'Passive Valving';


%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
simu.simMechanicsFile = 'Flap.slx';    % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer = 'on';                   % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = peakWavePeriod;                    % Wave Ramp Time [s]
simu.endTime = peakWavePeriod*20;                     % Simulation End Time [s]        
simu.solver = 'ode45';
simu.dt = 1e-1;                          % Simulation Time-Step [s]
simu.cicEndTime = 30;                   % Specify CI Time [s]
% simu.mcrMatFile = 'gains_mcr.mat';
% simu.mcrMatFile = 'pressure_mcr.mat';
% simu.mcrMatFile = 'PTO_wave_mcr.mat';


%% Wave Information
% % noWaveCIC, no waves with radiation CIC  
waves = waveClass('noWaveCIC');       % Initialize Wave Class and Specify Type  

% % Regular Waves 
% waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
% waves.height = 2.5;                     % Wave Height [m]
% waves.period = 8;                       % Wave Period [s]

% Irregular Waves using PM Spectrum with Directionality 
% waves = waveClass('irregular');         % Initialize Wave Class and Specify Type
if peakWavePeriod == 8  % wave height determined to set available wave power to 100kW (18*waves.power)
    waves.period = 8; waves.height = 1.22;
elseif peakWavePeriod == 20
    waves.period = 20; waves.height = 0.99;
end
waves.spectrumType = 'PM';              % Specify Spectrum Type
waves.phaseSeed = 1;
% waves.direction = [0,30,90];            % Wave Directionality [deg]
% waves.spread = [0.1,0.2,0.7];           % Wave Directional Spreading [%}

% % Irregular Waves with imported spectrum
% waves = waveClass('spectrumImport');      % Create the Wave Variable and Specify Type
% waves.spectrumFile = 'spectrumData.mat';  % Name of User-Defined Spectrum File [:,2] = [f, Sf]

% % Waves with imported wave elevation time-history  
% waves = waveClass('elevationImport');          % Create the Wave Variable and Specify Type
% waves.elevationFile = 'elevationData.mat';     % Name of User-Defined Time-Series File [:,2] = [time, eta]


%% Body Data
% Flap
body(1) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = 'geometry/flap.stl';     % Geometry File
body(1).mass = 127000;                          % User-Defined mass [kg]
body(1).inertia = [1.85e6 1.85e6 1.85e6];       % Moment of Inertia [kg-m^2]
body(1).setInitDisp([0 0 5], [0 1 0 .5], [0 0 0]); % Sets the initial rotation in y about the hinge rotation

% Base
body(2) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Base
body(2).geometryFile = 'geometry/base.stl';     % Geometry File
body(2).mass = 999;                             % Placeholder mass for a fixed body
body(2).inertia = [999 999 999];                % Placeholder inertia for a fixed body

%% PTO and Constraint Parameters
% Fixed
constraint(1)= constraintClass('Constraint1');  % Initialize ConstraintClass for Constraint1
constraint(1).location = [0 0 -10];             % Constraint Location [m]

% Rotational PTO
pto(1) = ptoClass('PTO1');                      % Initialize ptoClass for PTO1
pto(1).stiffness = 0;                           % PTO Stiffness Coeff [Nm/rad]
pto(1).damping = 0;%12000;                         % PTO Damping Coeff [Nsm/rad]
pto(1).location = [0 0 -8.9];                   % PTO Location [m]

%% load PTO parameters
params = Parameters();

%% Rotary to Linear Adjustable Rod
crank_z = 3; % [m] Seems to be the vertical distance up the flap that the rod is connected to
initalCylinderLength = 1.5*params.cylinderStroke;
crank_x = sqrt(initalCylinderLength^2 - crank_z^2);
% calculate the initial length of the cylinder
    % Fully retracted - the total length is about the stroke length
    % Fully extended - the total length is about 2x the stroke length
    % Half way - the total length is about 1.5x the stroke length
ptoSim(1) = ptoSimClass('adjustableRod');
ptoSim(1).adjustableRod.crank = crank_z;
ptoSim(1).adjustableRod.offset = 0;
ptoSim(1).adjustableRod.rodInit = 1.5*params.cylinderStroke; % [m] Initial length of the cylinder;
ptoSim(1).adjustableRod.initalCylinderLength = initalCylinderLength;


%% Parameters specific to case
params = UniqueParameters(PTO,peakWavePeriod,codesign,params);

function out = UniqueParameters(PTO,period,codesign,params)
% codesign only impacts the EHA PTO. It is ignored otherwise

% update params to reflect this specific case
out = params;

% Close to resonance
    closeToResonance = struct();
    % "Mechanical gains" reflect gains chosen by lookng at mechanicalenergy
        closeToResonance.mechGains.Kp = 5e5;
        closeToResonance.mechGains.Ki = -1e5;
    % "Electrical gains" reflect gains chosen by lookng at electricalenergy
        closeToResonance.elecGains.Kp = 6e6;
        closeToResonance.elecGains.Ki = 0;
    % Optimal pressure for Passive Valving PTO
        closeToResonance.rectifyingPressure = 5e6;

% Far from resonance
    farFromResonance = struct();
    % "Mechanical gains" reflect gains chosen by lookng at mechanicalenergy
        farFromResonance.mechGains.Kp = 4e6;
        farFromResonance.mechGains.Ki = 3e6;
    % "Electrical gains" reflect gains chosen by lookng at electricalenergy
        farFromResonance.elecGains.Kp = 3e6;
        farFromResonance.elecGains.Ki = 1e6;
    % Optimal pressure for Passive Valving PTO
        farFromResonance.rectifyingPressure = 16e6;

% Put these parameters into the params structure depending on the PTO chosen 
switch PTO
    case {'Active Valving', 'HHEA'}
        if period == 8
            out.Kp = farFromResonance.mechGains.Kp;
            out.Ki = farFromResonance.mechGains.Ki;
        elseif period == 20
            out.Kp = closeToResonance.mechGains.Kp;
            out.Ki = closeToResonance.mechGains.Ki;
        else
            disp('Gains not defined for this case')
            return
        end
    case 'EHA'
        if codesign
            if period == 8
                out.Kp = farFromResonance.elecGains.Kp;
                out.Ki = farFromResonance.elecGains.Ki;
            elseif period == 20
                out.Kp = closeToResonance.elecGains.Kp;
                out.Ki = closeToResonance.elecGains.Ki;
            else
                disp('Gains not defined for this case')
                return
            end
        else
            if period == 8
                out.Kp = farFromResonance.mechGains.Kp;
                out.Ki = farFromResonance.mechGains.Ki;
            elseif period == 20
                out.Kp = closeToResonance.mechGains.Kp;
                out.Ki = closeToResonance.mechGains.Ki;
            else
                disp('Gains not defined for this case')
                return
            end
        end % if codesign
    case 'Passive Valving'
        if period == 8
            out.pressure = farFromResonance.rectifyingPressure;
        elseif period == 20
            out.pressure = closeToResonance.rectifyingPressure;
        else
            disp('Gains not defined for this case')
            return
        end
    otherwise
        disp('PTO not defined')
        return
end % switch PTO
        
end % function

%% To do
    % Is this an okay wave condition? 100kW waves
    % Is the cylinder areas set okay so the pressure is reasonable?

    % Redo all the grid searches