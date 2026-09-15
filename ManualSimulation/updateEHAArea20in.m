% updateEHAArea20in.m
% One-off: updates results/sizedAreas.mat's EHA entry from the literal-
% peak-force sizing (23.37in bore, capArea=0.276822) to a 20in bore
% (capArea=0.202683 m^2, symmetric so rodArea=capArea too), per the
% peak-vs-RMS force discussion in diagnostics/ReadMe.md -- EHA's peak
% force (9.688MN) occurs only 0.093% of the time (RMS 2.595MN, peak/RMS
% ratio 3.73); a 20in bore (7.094MN capacity at 35MPa) clips the rare
% extreme spike but is much smaller than literal-peak sizing. Clipping
% itself is implemented in control/MPC_QP.m's EHA case (force saturation
% at capArea*hyd.maxPressure).
%
% Calls: none
% Called by: none (one-off top-level script)

repoRoot = fileparts(mfilename('fullpath'));
S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

oldCapArea = sizedAreas.EHA.capArea;
sizedAreas.EHA.capArea = 0.202683;
sizedAreas.EHA.rodArea = 0.202683;

save(fullfile(repoRoot,'results','sizedAreas.mat'),'sizedAreas');
fprintf('EHA capArea/rodArea updated: %.6f -> %.6f m^2 (20in bore)\n', oldCapArea, sizedAreas.EHA.capArea);
disp(sizedAreas);
