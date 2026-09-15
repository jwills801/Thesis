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
% same linearization already used throughout models/makeEHALossMap*.m).
% Like a gear-coupled flywheel, the shaft's own inertia J then reflects
% back to the flap side as J*k^2. Zero by default -- only applied when
% runParams explicitly selects the fixed-displacement EHA variant with an
% assumed shaft inertia (runParams.shaftInertia = 249.9 kgm^2 in every
% current call site: 28.6 kgm^2 from the hydraulic pump/motor side (one
% ~7297cc/rev pump, the actual displacement this sizing needs, trend-
% extrapolated -- J=a*D^b fit to real Rexroth A4VSO rotor-inertia-vs-
% displacement data -- rather than several smaller real units summed) +
% 221.3 kgm^2 from the electric motor/generator side (one 4000kW motor,
% sized so continuous rating covers roughly half of peak instantaneous
% power, with a 2x short-term overload assumed for the rest; also trend-
% extrapolated, from real ABB rotor-inertia-vs-power data), neither a
% real single catalog part -- see diagnostics/ReadMe.md's "EHA shaft
% inertia built up from real component data" section for the derivation).
phys.reflectedInertia = 0;
if nargin > 0 && isfield(runParams,'ehaFixedDisplacement') && runParams.ehaFixedDisplacement
    vMax = 1;
    w0 = 2000*(2*pi/60);
    Scale = vMax*runParams.capArea/w0*2*pi*1e6/107*1.2;
    d = (107*100^-3)/(2*pi);
    k = runParams.capArea*2.7574/(d*Scale);
    phys.reflectedInertia = runParams.shaftInertia * k^2;
end

% Linearized quadratic viscous damping (Lorentz/equivalent linearization).
% The TRUE nonlinear drag torque, from the Morison equation integrated
% over the flap's wetted height (Giorgi & Ringwood 2018, "Comparing
% nonlinear hydrodynamic forces in heaving point absorbers and
% oscillating wave surge converters", Eq.19: T_vis = sum_i L_i*(-0.5*rho*
% Cd*A_di*|V_i|*V_i), V_i=L_i*thetaDot), for a rectangular flap of width W
% rotating about a hinge at the seafloor, integrates in closed form to
%   T_vis = -(rho*Cd*W*Hwet^4/8) * |thetaDot|*thetaDot
% (a strip at height z has velocity z*thetaDot, frontal area W*dz, drag
% force -0.5*rho*Cd*W*z^2|thetaDot|thetaDot*dz, torque contribution z*dF;
% integrating z^3 from 0 to Hwet gives Hwet^4/4, times the -0.5*rho*Cd*W
% prefactor gives the Hwet^4/8 coefficient above). Hwet is the wetted
% height at rest -- the flap is 8.9m tall but the water is only 8m deep,
% so Hwet=8m (=phys.waterDepth), not the full flap height.
% Cd=4 is a dimensionless Morison drag coefficient -- checked against
% Giorgi & Ringwood 2018, who use Cd=8 for a comparable OWSC flap at low
% Keulegan-Carpenter number via Bearman et al. 1985's flat-plate formula;
% same order of magnitude, so this is a reasonable value (not
% independently re-derived here).
phys.Cd = 4; % dimensionless Morison drag coefficient
flapWidth = 18; % m
flapHeight = 8.9; % m
Hwet = min(flapHeight, phys.waterDepth); % m, wetted height at rest
bTrue = phys.rho * phys.Cd * flapWidth * Hwet^4 / 8; % Nm/(rad/s)^2 -- true nonlinear tau_v = bTrue*|thetaDot|*thetaDot coefficient
phys.bTrue = bTrue; % stored for later use by the true nonlinear plant (see note below); NOT the same as phys.b

% Lorentz-linearize about a representative velocity amplitude: for a
% sinusoidal thetaDot with amplitude V, the energy-equivalent LINEAR
% damping is B_lin=(8/(3*pi))*bTrue*V. V is taken as vMax/2.7574 -- the
% same vMax=1 m/s cylinder design-point velocity used throughout
% getHydraulic.m, converted to an angular velocity via the same
% cylinder-lever linearization (Force2Torque(0)=2.7574) used everywhere
% else in this codebase, as a representative "normal" operating
% amplitude. This is a modeling choice, not a measured value -- flag if
% a different nominal amplitude is wanted.
vMax = 1; % m/s
thetaDotNominal = vMax/2.7574; % rad/s
B_lin = (8/(3*pi)) * bTrue * thetaDotNominal;
% Baked directly into phys.sys.A below (not phys.b) so it applies
% identically to BOTH the controller's internal prediction (built from
% phys.sys via getTransition.m) and the true plant (advanceStep.m also
% propagates via sys.A) -- for now, both use this linearized damping;
% later the true plant can switch to the actual nonlinear
% bTrue*|thetaDot|*thetaDot term (via phys.b, still 0 by default) while
% the controller keeps the linear one. Setting phys.b=bTrue now would
% double-count damping in the true plant -- once via this linear
% A-matrix term, once via advanceStep.m's separate nonlinear
% viscousTorque term -- so phys.b is left at 0 here.

% State space
I_total = phys.I+phys.Iinf+phys.reflectedInertia;
A = [-B_lin/I_total, -phys.Khs/I_total, 0, -1/I_total;...
    1,               0,               0,        0;...
    0,               0,               0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0,               1, -x_timeDomain(2)];

B = [1/I_total;0;0;0]; 
C = eye(4);
phys.sys = ss(A,B,C,0);

end