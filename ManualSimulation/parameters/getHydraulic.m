function hyd = getHydraulic(runParams)

% Congifuration of hydraulic cylinder
hyd.rodArea = runParams.rodArea; % m^2
hyd.capArea = runParams.capArea; % m^2

% Cylinder torque and angular velocity
hyd.stroke = 5;
cylHorizDist = 7; % Distance from hinge to cylinder along sea floor
r_cyl = 3; % distance from hinge to cylinder along flap
hyd.L = @(theta) sqrt(cylHorizDist^2+r_cyl^2 + 2*cylHorizDist*r_cyl*sin(theta));
hyd.dLdt = @(theta,thetaDot) cylHorizDist*r_cyl*cos(theta).*thetaDot./hyd.L(theta);
hyd.Force2Torque = @(theta) cylHorizDist*r_cyl*cos(theta)./hyd.L(theta);

% Find Important Lengths
hyd.L_equilib = hyd.L(0);
hyd.L_retract = hyd.L_equilib - hyd.stroke/2; % Length of cylinder at full retraction

switch runParams.drive
    case {'DHD', 'PassivePump'}
        % Set Valued Control Inputs
        hyd.pressureRails = runParams.pressureRails;
        capForceOptions = hyd.pressureRails * hyd.capArea;
        rodForceOptions = hyd.pressureRails * hyd.rodArea;
        hyd.ptoForceOptions = capForceOptions'-rodForceOptions;
end


% load or calculate switching losses
switch runParams.drive
    case 'DHD'
        switchMap = makeSwitchLossMap(hyd);
        save("parameters/SwitchMap.mat","switchMap")
        % load("SwitchMap.mat")
        hyd.switchMap = switchMap;
    case 'EHA'
        % Load EHA losses
        vMax = 1;
        hyd.EHA = makeEHALossMap(hyd.capArea,vMax);
end

% Output function handels
hyd.getVolandFlow = @getVolandFlow;

end

% Functions
function [cap,rod] = getVolandFlow(params,states)
    theta = states(2,:)';
    thetaDot = states(1,:)';
    
    % Calculate the length of the whole cylinder
    L = params.hyd.L(theta);
    
    % Calcuate the disctance from TDC
    xCap = L-params.hyd.L_retract;
        
    % Calulcate the distance from BDC
    xRod = params.hyd.stroke - xCap;

    % Calculate volumes in each side
    cap.vol = xCap*params.hyd.capArea;
    rod.vol = xRod*params.hyd.rodArea;

    % Calculate velocity of the cylinder
        % positive is in extension
    dLdt = params.hyd.dLdt(theta,thetaDot);

    % Calculate ideal flow into each cylinder
    cap.velA = dLdt*params.hyd.capArea;
    rod.velA = -dLdt*params.hyd.rodArea;
end



