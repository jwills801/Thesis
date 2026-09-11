% getSimulation.m
% Default simulation-timing/sea-state parameters (finalTime, dt,
% peakPeriod, sigWaveHeight, makePlots flag). Ignores its argument;
% callers override individual fields on the returned struct afterward
% (e.g. params.simu.sigWaveHeight = seaState.Hs) rather than passing them
% in here.
% Calls: none
% Called by: parameters/getParameters.m
function simu = getSimulation(~)

simu.finalTime = 500;
simu.waveFinalTime = 600;
simu.dt = 1e-2;
simu.rampTime = 50; % s
simu.time = (0:simu.dt:simu.finalTime)';

simu.peakPeriod = 8; % s
simu.sigWaveHeight = 2.5; % [m]

% simu.waveType = "monochromatic";
simu.waveType = "polychromatic";

% Set to false for batch/sweep runs (e.g. diagnostics, Run_All_Cases) to
% suppress the diagnostic figures some evaluation functions pop up.
simu.makePlots = true;

end