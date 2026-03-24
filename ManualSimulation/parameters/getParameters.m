function params = getParameters(~)
params = struct;
params.phys = getPhysical;
params.hyd = getHydraulic;
params.simu = getSimulation;
end