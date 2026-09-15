% Run_EHAHydraulicElectricEfficiencyMaps.m
% Decomposes the EHA's fixed-displacement pump-motor loss physics
% (models/makeEHALossMap_fixedDisp.m's ehaLoss, the REAL physics -- not
% the least-squares quadratic fit EHA.LossCoeffs) into its two stages and
% plots the EFFICIENCY of each, side by side, using the explicit
% motoring/generating branch agreed on (not the min/max shortcut):
%
%   P_hyd = Q*deltaP; P_shaft = w*T_Act; P_elec = P_shaft + P_L_elect
%   (P_L_elect = copperCoeff*T_Act^2, unsigned -- the sign(T_Act) that
%   used to be here was a mistaken direction proxy, removed from
%   models/makeEHALossMap_fixedDisp.m; see diagnostics/ReadMe.md)
%   Direction: P_hyd>0 -> motoring (grid->elec->shaft->hyd->flap)
%              P_hyd<0 -> generating (flap->hyd->shaft->elec->grid)
%   Motoring:   hyd_eff = |P_hyd|/|P_shaft|,   elec_eff = |P_shaft|/|P_elec|
%   Generating: hyd_eff = |P_shaft|/|P_hyd|,   elec_eff = |P_elec|/|P_shaft|
%
% Deliberately NOT clamped to <=1 -- efficiency >1 at a grid point is a
% visible flag that the model's output magnitude exceeded its input
% there (the known sign-convention inconsistency), not hidden by a
% min/max clamp. Same grid as plotting/plotEHAEfficiencyMap.m.
%
% Calls: none
% Called by: none (top-level entry point)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'parameters'), fullfile(repoRoot,'models'));

S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
capArea = S.sizedAreas.EHA.capArea;
vMax = 1;
copperCoeff = 5.0e-4;

% Same grid as plotting/plotEHAEfficiencyMap.m
nW = 60; W_vals = linspace(-.5,.5,nW);
nT = 60; T_vals = linspace(-1,1,nT)*3e7;
[Wg,Tg] = ndgrid(W_vals,T_vals);

hydEff = NaN(size(Wg));
elecEff = NaN(size(Wg));

Wrpm = 2000; w0 = Wrpm*(2*pi/60);
Scale = vMax*capArea/w0*2*pi*1e6/107*1.2;
D = 107; d = (D*100^-3)/(2*pi);
Cf =  53.7e-3; Ch = 53.6; Cv = 23.5e3; Cs = 4.26e-9;
mu=(32e-6)*870; B = 1.7e9; rho = 870;

for i = 1:numel(Wg)
    thetaDot = Wg(i); T_flap = Tg(i);
    V = thetaDot*2.7574; Q = V*capArea;
    F = T_flap/2.7574; deltaP = F/capArea;

    QLoss = Scale*abs(d*Cs*deltaP/mu) + Scale*abs(d*w0*deltaP/B);
    w = (Q + sign(deltaP)*QLoss) / (d*Scale);

    TLoss = Scale*(abs(d*Cv*mu*w) + abs(d*deltaP*Cf) + abs(Ch*w^2*rho*d^(5/3)/2));
    T_Ideal = deltaP*d*Scale;
    T_Act = T_Ideal + sign(w)*TLoss;

    P_hyd = Q*deltaP;
    P_shaft = w*T_Act;
    P_L_elect = copperCoeff*T_Act^2;
    P_elec = P_shaft + P_L_elect;

    if P_hyd > 0 % motoring: grid -> elec -> shaft -> hyd -> flap
        hydEff(i) = abs(P_hyd)/max(abs(P_shaft),1);
        elecEff(i) = abs(P_shaft)/max(abs(P_elec),1);
    else % generating (or exactly zero): flap -> hyd -> shaft -> elec -> grid
        hydEff(i) = abs(P_shaft)/max(abs(P_hyd),1);
        elecEff(i) = abs(P_elec)/max(abs(P_shaft),1);
    end
end

fprintf('hydEff range (unclipped): [%.3f, %.3f]  (values >1 flag the known sign inconsistency)\n', min(hydEff(:)), max(hydEff(:)));
fprintf('elecEff range (unclipped): [%.3f, %.3f]\n', min(elecEff(:)), max(elecEff(:)));
fprintf('%d/%d hydEff points and %d/%d elecEff points exceed 1 (clipped for display below)\n', ...
    sum(hydEff(:)>1), numel(hydEff), sum(elecEff(:)>1), numel(elecEff));

%% Plot -- one figure, side by side, same style as plotEHAEfficiencyMap.m's Figure 2
% Display clipped to [0,1] (data itself is untouched above/in the printed
% ranges) -- values >1 are a known sign-convention artifact, not a
% meaningful super-unity efficiency, so showing them isn't informative.
hydLevels = [.2 .5 .7 .8 .85 .9 .95 1];
elecLevels = linspace(0,1,11);
hydEffDisp = min(hydEff,1);
elecEffDisp = min(elecEff,1);

fig = figure('Position',[100 100 1300 550],'Visible','off');
subplot(1,2,1)
contourf(Wg,Tg/1e6,hydEffDisp,hydLevels,'ShowText','on');
colorbar, xlabel('\theta_{dot} [rad/s]'), ylabel('Torque [MNm]')
title('Hydraulic pump/motor efficiency')

subplot(1,2,2)
contourf(Wg,Tg/1e6,elecEffDisp,elecLevels,'ShowText','on');
colorbar, xlabel('\theta_{dot} [rad/s]'), ylabel('Torque [MNm]')
title('Electric machine efficiency')

sgtitle('EHA efficiency by stage, explicit motoring/generating branch (copperCoeff=5.0e-4)');
print(fig, fullfile(repoRoot,'results','ehaStageEfficiencyMaps.png'), '-dpng', '-r150', '-painters');
fprintf('Saved.\n');
