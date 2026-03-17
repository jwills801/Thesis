function plotting(params,wave,optTraj,dyn,eval)

% Plot energy loss per event
figure, yyaxis right
plot(eval.switchTimes(1:end-1),eval.lossAtSwitches,'*'), ylabel('Valve Loss [W]')
yyaxis left
plot(dyn.t,dyn.u), xlabel('Time [s]'), ylabel('Control Input [Nm]'), grid

% Plot cummulative valve loss over time
figure, plot(eval.switchTimes(1:end-1),cumsum(eval.lossAtSwitches),...
    dyn.t,eval.mechEnergy), legend('Switching Loss','Absorbed Energy')
xlabel('Time [s]'), ylabel('Cummulative Energy [J]'), grid

% Plot control performance
figure
subplot(221), plot(dyn.t,dyn.u), xlabel('Time [s]'), ylabel('Control Input [Nm]'), grid
subplot(222), yyaxis left, plot(dyn.t,dyn.thetaDot), ylabel('Angular Velocity [rad/s]'), grid
yyaxis right,plot(wave.torque.time,wave.torque.Texc), xlabel('Time [s]'), ylabel('Excitaiton Torque [Nm]')
subplot(223), plot(dyn.t,dyn.thetaDot,'k',optTraj.time,optTraj.thetaDot,'k--'), xlabel('Time [s]'), ylabel('Angular Velocity [W]'), legend('Actual','Optimal'), grid
% subplot(224), plot(t,[0,s]), xlabel('Time [s]'), ylabel('Sliding Surface'), grid

% Compare average powers
disp(['Optimal Power = ', num2str(optTraj.avePow/1e3,3), 'kW']);
disp(['Mechanical Power = ', num2str(eval.aveMechPow/1e3,3), ' kW']);
disp(['Electric Power = ', num2str(eval.aveElecPow/1e3,3), ' kW']);
end