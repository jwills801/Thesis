% getParameters.m
% Assembles the top-level params struct (phys/hyd/simu/runParams) for a
% run from runParams.
% Calls: getPhysical.m, getHydraulic.m, getSimulation.m
% Called by: main_WEC_Simulation.m, parameters/optimizePressure.m,
%   sizeCylinderArea.m, diagnostics/checkAstarVsBruteForce.m,
%   checkEnergyBalance.m, checkMPC_EHA.m, checkTerminalCost.m,
%   checkTransitionMatrices.m, validatePhase2Subset.m
function params = getParameters(runParams)
params = struct;
params.phys = getPhysical(runParams);
params.hyd = getHydraulic(runParams);
params.simu = getSimulation;
params.runParams = runParams;
end