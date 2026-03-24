function simu = getSimulation(~)

simu.finalTime = 100;
simu.dt = 1e-2;
simu.rampTime = 50; % s
simu.time = (0:simu.dt:simu.finalTime)';

simu.peakPeriod = 8; % s
simu.sigWaveHeight = 2.5; % [m]

end