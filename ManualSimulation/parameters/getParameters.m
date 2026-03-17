function params = getParameters(~)
params = struct;
params.phys = getPhysical;
params.hyd = getHydraulic;
params.cntrl = getControl;
params.simu = getSimulation;
end