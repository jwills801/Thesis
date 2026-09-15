% checkPhysicalStateSpace.m
% Rebuilds A,B independently from the main(1).tex "Full State Space"
% derivation and diffs against params.phys.sys, to catch transcription
% bugs in getPhysical.m. Calls getPhysical() with no runParams, so it
% only checks the baseline (zero reflected-inertia) case.
% Calls: parameters/getPhysical.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'parameters'));

phys = getPhysical();

% Rebuild from the LaTeX derivation directly (independent arithmetic from
% getPhysical.m's own code, same source numbers)
rho = 1024; Vol = 297; g = 9.81; r_cob = 4; r_cog = 5; m = 127000;
Khs = (rho*Vol*r_cob - m*r_cog)*g;                 % K_hs = (rho*V*r_cob - m*r_cog)*g
I = 5.025e6; Iinf = 1.734e7;
x1 = 1.9754; x2 = 1.1345; x3 = 7.6921;              % x*_t from main(1).tex
Itot = I + Iinf;

% Lorentz-linearized quadratic damping: B_lin = (8/(3*pi))*bTrue*V, with
% bTrue the Morison-integrated true nonlinear coefficient (rho*Cd*W*
% Hwet^4/8) and V=vMax/2.7574 (vMax=1 m/s, same design point used in
% getHydraulic.m). See getPhysical.m's comment for the full derivation.
rhoW = 1024; Cd = 4; flapWidth = 18; flapHeight = 8.9; waterDepthLocal = 8;
Hwet = min(flapHeight, waterDepthLocal);
bTrue = rhoW*Cd*flapWidth*Hwet^4/8;
vMax = 1; thetaDotNominal = vMax/2.7574;
Blin = (8/(3*pi))*bTrue*thetaDotNominal;

A_expected = [-Blin/Itot, -Khs/Itot, 0, -1/Itot;
              1,         0,        0,  0;
              0,         0,        0, -x1;
              x3*1e7,    0,        1, -x2];
B_expected = [1/Itot; 0; 0; 0];

dKhs = abs(Khs - phys.Khs);
dA = max(abs(A_expected(:) - phys.sys.A(:)));
dB = max(abs(B_expected(:) - phys.sys.B(:)));

fprintf('--- checkPhysicalStateSpace ---\n');
report('Khs matches (rho*V*r_cob - m*r_cog)*g', dKhs, 1e-6);
report('A matrix matches Full State Space derivation', dA, 1e-9);
report('B matrix matches Full State Space derivation', dB, 1e-9);

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g)\n',name,err);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
