function simu = getSimulation(~)

simu.finalTime = 50;
simu.dt = 1e-2;
simu.rampTime = 1e-3; % s
simu.time = (0:simu.dt:simu.finalTime)';

simu.peakPeriod = 10; % s
simu.sigWaveHeight = 4; % [m]

end