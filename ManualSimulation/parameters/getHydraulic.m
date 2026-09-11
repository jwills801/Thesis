% getHydraulic.m
% Builds hyd: cylinder geometry, pressure rails (evenlySpacedRails.m for
% 3/4-rail DHD/PassivePump interior placement), PTO force options, the
% DHD switch-loss map (makeSwitchLossMap.m, or PassivePump's check-valve
% constant), and the EHA loss map (makeEHALossMap.m, or
% makeEHALossMap_fixedDisp.m if runParams.ehaFixedDisplacement).
% Calls: evenlySpacedRails.m, makeSwitchLossMap.m, makeEHALossMap.m,
%   makeEHALossMap_fixedDisp.m, plotEHAEfficiencyMap.m (opt-in diagnostic
%   plot only, gated behind runParams.plotEHAEfficiency)
% Called by: parameters/getParameters.m
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

% Reservoir pressure -- always 0.5MPa regardless of highPressure (fixing
% a previous mistake where the low rail was scaled proportionally with
% highPressure, e.g. giving 0.29MPa at highPressure=20.6MPa instead of
% the physically-fixed reservoir pressure).
p0 = .5; % MPa, fixed reservoir pressure
switch runParams.drive
    case {'DHD', 'PassivePump'}
        pLow = p0*1e6;
        pHigh = runParams.highPressure;
        switch runParams.pressure_rails
            case 2
                hyd.pressureRails = [pLow, pHigh];
            case 3
                % Interior (middle) rail placed to make the resulting DHD
                % force options as evenly spaced as possible, rather than
                % the previous fixed fraction of highPressure.
                hyd.pressureRails = evenlySpacedRails(pLow,pHigh,hyd.capArea,hyd.rodArea,3);
            case 4
                % Both interior rails placed the same way.
                hyd.pressureRails = evenlySpacedRails(pLow,pHigh,hyd.capArea,hyd.rodArea,4);
            case 5
                % Unchanged fixed-fraction interior placement (out of the
                % current comparison's scope, only 2/3/4 rails are being
                % compared) but same low-rail fix: interior/top rails keep
                % their original relative position between pLow and pHigh
                % instead of also being scaled down with highPressure.
                oldRails = [3.5 18.5 29 35]*1e6;
                frac = (oldRails - p0*1e6) / (35e6 - p0*1e6);
                hyd.pressureRails = [pLow, pLow + frac*(pHigh-pLow)];
        end

        capForceOptions = hyd.pressureRails * hyd.capArea;
        rodForceOptions = hyd.pressureRails * hyd.rodArea;
        hyd.ptoForceOptions = capForceOptions'-rodForceOptions;
end


% Design-point cylinder velocity, used both for EHA loss-map sizing below
% and for the PassivePump check-valve sizing.
vMax = 1; % m/s

% load or calculate switching losses
switch runParams.drive
    case 'DHD'
        %switchMap = makeSwitchLossMap(hyd); save("parameters/SwitchMap.mat","switchMap")
        load("SwitchMap.mat")
        hyd.switchMap = switchMap;
    case 'PassivePump'
        % getValveLoss.m's PassivePump branch needs switchMap.valveConstant
        % for its steady open-valve loss (this used to reuse the DHD
        % switching valve's spec -- "2WRC-4x size 80" -- via a shared
        % SwitchMap.mat, which both crashed when that case narrowed to
        % 'DHD'-only, and physically didn't make sense: PassivePump has a
        % fixed check valve, not a fast-switching valve, so it shouldn't
        % share that sizing). Size the check valve from the orifice
        % equation Q=valveConstant*sqrt(dP) so it drops ~0.5MPa at the
        % worst-case (larger, cap-side) max flow, using the same vMax=1
        % m/s cylinder-velocity design point used for EHA sizing below.
        checkValveDP = 0.5e6; % Pa
        Qmax = vMax*hyd.capArea;
        hyd.switchMap.valveConstant = Qmax/sqrt(checkValveDP);
end
% Load EHA losses. Fixed speed (2000 RPM, variable displacement) is the
% default/baseline; set runParams.ehaFixedDisplacement=true (plus
% runParams.shaftInertia, consumed by getPhysical.m for the reflected
% inertia) to switch to fixed displacement (chi=1), variable speed.
if isfield(runParams,'ehaFixedDisplacement') && runParams.ehaFixedDisplacement
    hyd.EHA = makeEHALossMap_fixedDisp(hyd.capArea,vMax);
else
    hyd.EHA = makeEHALossMap(hyd.capArea,vMax);
end

% Opt-in diagnostic: set runParams.plotEHAEfficiency=true to visualize the
% EHA loss map (real pump physics vs. the quadratic fit the electrical-
% optimized MPC_QP controller uses internally). Off by default -- this is
% a diagnostic plot, not something every simulation run should pop up.
if isfield(runParams,'plotEHAEfficiency') && runParams.plotEHAEfficiency
    plotEHAEfficiencyMap(hyd.EHA, hyd.capArea);
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



