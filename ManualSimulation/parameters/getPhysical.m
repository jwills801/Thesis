% getPhysical.m
% Builds the flap's physical parameters and 4-state linear state-space
% model (phys.sys): hydrostatic stiffness, radiation/added inertia,
% quadratic viscous damping coefficient (phys.b, zero by default), and
% the reflected inertia from a fixed-displacement EHA (phys.reflectedInertia,
% zero unless runParams.ehaFixedDisplacement is set).
% Calls: none
% Called by: parameters/getParameters.m; diagnostics/checkPhysicalStateSpace.m
%   calls it directly (with no runParams) to check the baseline (zero
%   reflected-inertia) case.
function phys = getPhysical(runParams)


phys.rho = 1024;
phys.Vol = 297;
phys.g = 9.81;
phys.r_cob = 4;
phys.r_cog = 5;
phys.waterDepth = 8;
phys.m = 127000;

% Quadratic viscous damping coefficient b in tau_v(t) = b|thetaDot|thetaDot
% (main(1).tex "Inclusion of Viscous Damping"). Zero by default: none of
% the existing linear-model results include any viscous damping term at
% all. dynamics/advanceStep.m applies this to the TRUE plant only --
% controllers' internal prediction models (getControl.m, getOptimal.m)
% stay built on the zero-damping linear system regardless of this value.
% Set explicitly (e.g. via runParams) for the phase-3 nonlinear-damping
% robustness sweep; there's no existing value in this codebase to anchor
% a nonzero default to.
phys.b = 0;

% hydrostatic stiffness
phys.Khs = phys.rho * phys.Vol * phys.g * phys.r_cob - ...
    phys.m * phys.g * phys.r_cog;

% Inertia
phys.I = 5.025e6;

% Radiation
phys.Iinf = 1.734e7;
x_timeDomain = [1.9754, 1.1345, 7.6921];

% Reflected inertia from a fixed-displacement (chi=1) EHA. Flow continuity
% between cylinder and pump forces shaft speed omega = k*thetaDot, with
% k = capArea*2.7574/(d*Scale) the pump-to-shaft "gear ratio" (2.7574 is
% Force2Torque(0), the cylinder-lever ratio linearized at theta=0, the
% same linearization already used throughout parameters/makeEHALossMap*.m).
% Like a gear-coupled flywheel, the shaft's own inertia J then reflects
% back to the flap side as J*k^2. Zero by default -- only applied when
% runParams explicitly selects the fixed-displacement EHA variant with an
% assumed shaft inertia (runParams.shaftInertia; no measured value exists
% for this system, see diagnostics/ReadMe.md).
phys.reflectedInertia = 0;
if nargin > 0 && isfield(runParams,'ehaFixedDisplacement') && runParams.ehaFixedDisplacement
    vMax = 1;
    w0 = 2000*(2*pi/60);
    Scale = vMax*runParams.capArea/w0*2*pi*1e6/107*1.2;
    d = (107*100^-3)/(2*pi);
    k = runParams.capArea*2.7574/(d*Scale);
    phys.reflectedInertia = runParams.shaftInertia * k^2;
end

% State space
I_total = phys.I+phys.Iinf+phys.reflectedInertia;
A = [0,            -phys.Khs/I_total, 0, -1/I_total;...
    1,               0,               0,        0;...
    0,               0,               0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0,               1, -x_timeDomain(2)];

B = [1/I_total;0;0;0]; 
C = eye(4);
phys.sys = ss(A,B,C,0);

end