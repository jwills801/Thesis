function params = getParameters(runParams)
params = struct;
params.phys = getPhysical;
params.hyd = getHydraulic(runParams);
params.simu = getSimulation;
params.runParams = runParams;
end