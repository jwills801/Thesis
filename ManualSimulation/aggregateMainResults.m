function aggregateMainResults()
% aggregateMainResults.m
% Standalone aggregator for Run_MainResults.m -- deliberately NOT called
% automatically at the end of that script, so it can be run at any time
% against whatever results/mainResults/task_*.mat files exist on disk,
% including while the sweep is still running or stuck on a few hung
% pressure points (see Run_MainResults.m's header for why tasks are split
% per-pressure-point rather than per-family/sea-state).
%
% For PassivePump/DHD2/DHD3/DHD4, rebuilds each (family, sea state)'s
% INTENDED pressure grid (same prior-lookup + fine/coarse logic as
% Run_MainResults.m) so coverage can be reported as "found/intended", not
% just a raw count. Selection, applied to EVERY rail family (not just
% PassivePump): argmax(elecRGP) among finished points that swing
% SYMMETRICALLY past a +/-5deg position target on both sides
% (maxThetaDeg>5 AND minThetaDeg<-5 -- not just peak|theta| one-sided,
% which a lopsided oscillation could satisfy without the flap actually
% swinging through a healthy range in both directions; this matters
% because of standing concerns about Coulomb-damping simulation validity
% at large damping magnitudes), falling back to whichever already-run
% pressure is closest to 5MPa (the floor -- never a new simulation) if no
% pressure in that combo's grid clears the symmetric bar. PassivePump
% SS1/SS3 needed a real 5MPa point added once (see
% Run_PassivePumpSS3_5MPa.m) since their fine grids didn't otherwise
% reach it; DHD families use whatever's already closest on disk instead
% of requiring an exact 5MPa point, per instruction not to run new sims
% for this.
%
% Writes results/mainResults.csv.
%
% Calls: none
% Called by: none (standalone, rerun anytime)

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'wave'));
resultsDir = fullfile(repoRoot,'results','mainResults');

seaStates = humboldtSeaStates();
nSeaStates = height(seaStates);
POSITION_TARGET_DEG = 5;

FINE_OFFSETS = [-11 -6 -4 -2 0 2 4 6 11] * 1e6;
COARSE_GRID = [5 11 17 23 29 35] * 1e6;

priorFile = fullfile(repoRoot,'results','mAstarConvAndPressureOpt','pressureOptSummary.csv');
priorTable = readtable(priorFile);
priors = containers.Map('KeyType','char','ValueType','double');
for i = 1:height(priorTable)
    key = sprintf('%s_%d', priorTable.Family{i}, priorTable.SeaStateIdx(i));
    priors(key) = priorTable.BestPressure_MPa(i) * 1e6;
end
priors('DHD4_4') = 20e6;

railFamilies = {'PassivePump','DHD2','DHD3','DHD4'};

rows = struct('Family',{},'SeaStateIdx',{},'Hs',{},'Tp',{},'Probability',{}, ...
    'BestPressure_MPa',{},'MechRGP',{},'ElecRGP',{},'AveMechPow_kW',{},'AveElecPow_kW',{}, ...
    'MaxThetaDeg',{},'MinThetaDeg',{},'MetPositionTarget',{},'PointsFound',{},'PointsIntended',{});

for f = 1:numel(railFamilies)
    family = railFamilies{f};
    for ss = 1:nSeaStates
        Hs = seaStates.Hs(ss); Tp = seaStates.Tp(ss); prob = seaStates.probability(ss);
        key = sprintf('%s_%d', family, ss);
        if isKey(priors,key)
            grid = unique(min(max(priors(key) + FINE_OFFSETS, 5e6), 35e6));
        else
            grid = COARSE_GRID;
        end
        if strcmp(family,'PassivePump')
            % Always include the 5MPa floor as a candidate, even if the
            % fine grid's offsets don't reach it -- it's the designated
            % fallback point when no pressure clears the symmetric
            % +/-10deg position target (see the PassivePump selection
            % logic below), so it must always be checked for/counted.
            grid = unique([grid, 5e6]);
        end
        nIntended = numel(grid);

        found = struct('pressure',{},'mechRGP',{},'elecRGP',{},'aveMechPow',{},'aveElecPow',{}, ...
            'maxThetaDeg',{},'minThetaDeg',{});
        for p = grid
            pMPa = p/1e6;
            fname = fullfile(resultsDir, sprintf('task_%s_ss%d_p%05.1f.mat', family, ss, pMPa));
            if isfile(fname)
                L = load(fname);
                found(end+1) = struct('pressure',L.out.pressure,'mechRGP',L.out.mechRGP, ...
                    'elecRGP',L.out.elecRGP,'aveMechPow',L.out.aveMechPow,'aveElecPow',L.out.aveElecPow, ...
                    'maxThetaDeg',L.out.maxThetaDeg,'minThetaDeg',L.out.minThetaDeg); %#ok<AGROW>
            end
        end
        nFound = numel(found);

        if nFound == 0
            rows(end+1) = struct('Family',family,'SeaStateIdx',ss,'Hs',Hs,'Tp',Tp,'Probability',prob, ...
                'BestPressure_MPa',NaN,'MechRGP',NaN,'ElecRGP',NaN,'AveMechPow_kW',NaN,'AveElecPow_kW',NaN, ...
                'MaxThetaDeg',NaN,'MinThetaDeg',NaN,'MetPositionTarget',NaN, ...
                'PointsFound',0,'PointsIntended',nIntended); %#ok<AGROW>
            continue
        end

        elecRGPs = [found.elecRGP];
        maxThetas = [found.maxThetaDeg];
        minThetas = [found.minThetaDeg];
        pressuresFound = [found.pressure];

        % Symmetric +/-5deg position target, applied to every rail family
        % (not just PassivePump anymore): among points that clear it on
        % BOTH sides, pick argmax(elecRGP); otherwise fall back to
        % whichever already-run pressure is closest to 5MPa (the floor)
        % -- never a new simulation, just the closest point already on
        % disk for that combo.
        candidates = find(maxThetas > POSITION_TARGET_DEG & minThetas < -POSITION_TARGET_DEG);
        if isempty(candidates)
            [~,bestInd] = min(abs(pressuresFound - 5e6));
        else
            [~,bestLocal] = max(elecRGPs(candidates));
            bestInd = candidates(bestLocal);
        end
        metTarget = double(maxThetas(bestInd) > POSITION_TARGET_DEG && minThetas(bestInd) < -POSITION_TARGET_DEG);

        b = found(bestInd);
        rows(end+1) = struct('Family',family,'SeaStateIdx',ss,'Hs',Hs,'Tp',Tp,'Probability',prob, ...
            'BestPressure_MPa',b.pressure/1e6,'MechRGP',b.mechRGP,'ElecRGP',b.elecRGP, ...
            'AveMechPow_kW',b.aveMechPow/1e3,'AveElecPow_kW',b.aveElecPow/1e3, ...
            'MaxThetaDeg',b.maxThetaDeg,'MinThetaDeg',b.minThetaDeg,'MetPositionTarget',metTarget, ...
            'PointsFound',nFound,'PointsIntended',nIntended); %#ok<AGROW>
    end
end

% EHA fixed-displacement, electric-loss-aware -- no grid
for ss = 1:nSeaStates
    Hs = seaStates.Hs(ss); Tp = seaStates.Tp(ss); prob = seaStates.probability(ss);
    fname = fullfile(resultsDir, sprintf('task_EHA_fixed_elec_ss%d.mat', ss));
    if isfile(fname)
        L = load(fname);
        rows(end+1) = struct('Family','EHA_fixed_elec','SeaStateIdx',ss,'Hs',Hs,'Tp',Tp,'Probability',prob, ...
            'BestPressure_MPa',NaN,'MechRGP',L.out.mechRGP,'ElecRGP',L.out.elecRGP, ...
            'AveMechPow_kW',L.out.aveMechPow/1e3,'AveElecPow_kW',L.out.aveElecPow/1e3, ...
            'MaxThetaDeg',NaN,'MinThetaDeg',NaN,'MetPositionTarget',NaN, ...
            'PointsFound',1,'PointsIntended',1); %#ok<AGROW>
    else
        rows(end+1) = struct('Family','EHA_fixed_elec','SeaStateIdx',ss,'Hs',Hs,'Tp',Tp,'Probability',prob, ...
            'BestPressure_MPa',NaN,'MechRGP',NaN,'ElecRGP',NaN,'AveMechPow_kW',NaN,'AveElecPow_kW',NaN, ...
            'MaxThetaDeg',NaN,'MinThetaDeg',NaN,'MetPositionTarget',NaN, ...
            'PointsFound',0,'PointsIntended',1); %#ok<AGROW>
    end
end

T = struct2table(rows);
outFile = fullfile(repoRoot,'results','mainResults.csv');
writetable(T, outFile);
fprintf('Wrote %s\n', outFile);

%% Coverage + probability-weighted summary
fprintf('\n--- Coverage (found/intended pressure points, sea states with data) ---\n');
allFamilies = [railFamilies, {'EHA_fixed_elec'}];
for f = 1:numel(allFamilies)
    family = allFamilies{f};
    fam = T(strcmp(T.Family,family),:);
    ssWithData = sum(fam.PointsFound > 0);
    totalFound = sum(fam.PointsFound); totalIntended = sum(fam.PointsIntended);
    fprintf('%-16s %d/%d sea states have data, %d/%d pressure points total\n', ...
        family, ssWithData, nSeaStates, totalFound, totalIntended);
end

fprintf('\n--- Probability-weighted (expected, across the wave climate; NaN sea states dropped, not zero-filled) ---\n');
for f = 1:numel(allFamilies)
    family = allFamilies{f};
    fam = T(strcmp(T.Family,family) & ~isnan(T.ElecRGP),:);
    if isempty(fam), fprintf('%-16s no data yet\n', family); continue; end
    normProb = fam.Probability / sum(fam.Probability); % renormalize over available sea states
    expMech = sum(normProb .* fam.MechRGP);
    expElec = sum(normProb .* fam.ElecRGP);
    fprintf('%-16s expected mechRGP=%.4f expected elecRGP=%.4f (over %d/%d sea states with data)\n', ...
        family, expMech, expElec, height(fam), nSeaStates);
end
end
