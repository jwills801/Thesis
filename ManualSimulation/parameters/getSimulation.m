function simu = getSimulation(~)

simu.finalTime = 16;
simu.waveFinalTime = 50;
simu.dt = 1e-2;
simu.rampTime = 8; % s
simu.time = (0:simu.dt:simu.finalTime)';

simu.peakPeriod = 8; % s
simu.sigWaveHeight = 2.5; % [m]

end