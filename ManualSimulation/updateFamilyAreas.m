% updateFamilyAreas.m
% One-off: locks in the bore diameters chosen from the cylinder-area
% sweep (results/cylinderAreaSweep/summary.csv, sea state 8, 35MPa,
% mAstar=5 for DHD): PassivePump=13in, DHD2=16in, DHD3=17in, DHD4=18in.
% Each family now has its OWN distinct bore (previously DHD2/3/4 shared
% one area) -- capArea/rodArea use the same 1.5:1 ratio convention as
% every other sizing script in this repo. EHA's area (20in, sized
% separately) is left untouched.
% Calls: none
% Called by: none (one-off, already run)
repoRoot = fileparts(mfilename('fullpath'));
S = load(fullfile(repoRoot,'results','sizedAreas.mat'));
sizedAreas = S.sizedAreas;

in2m = 0.0254;
diamToAreas = @(d_in) deal(pi*(d_in*in2m/2)^2, pi*(d_in*in2m/2)^2/1.5);

diams = struct('PassivePump',13, 'DHD2',16, 'DHD3',17, 'DHD4',18);
fn = fieldnames(diams);
for i = 1:numel(fn)
    fam = fn{i};
    [capArea, rodArea] = diamToAreas(diams.(fam));
    sizedAreas.(fam).capArea = capArea;
    sizedAreas.(fam).rodArea = rodArea;
    fprintf('%-12s %2gin -> capArea=%.6f rodArea=%.6f\n', fam, diams.(fam), capArea, rodArea);
end

save(fullfile(repoRoot,'results','sizedAreas.mat'),'sizedAreas');
fprintf('Saved.\n');
