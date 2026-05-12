function eval = plotAll(params,wave,ctrl,dyn,eval)



% Plot control performance
figure
subplot(221), plot(dyn.t,dyn.u), xlabel('Time [s]'), ylabel('Control Input [Nm]'), grid
subplot(222), yyaxis left, plot(dyn.t,dyn.thetaDot), ylabel('Angular Velocity [rad/s]'), grid
yyaxis right,plot(wave.torque.time,wave.torque.Texc), xlabel('Time [s]'), ylabel('Excitaiton Torque [Nm]'), grid, xlim([min(dyn.t) max(dyn.t)])
subplot(223), plot(dyn.t,dyn.thetaDot,'k',ctrl.optTraj.time,ctrl.optTraj.thetaDot,'k--'), xlabel('Time [s]'), ylabel('Angular Velocity [W]'), legend('Actual','Optimal'), grid, , xlim([min(dyn.t) max(dyn.t)])
switch params.runParams.controller
    case 'PI'
        
    case 'Sliding Mode'
        thetaErr = dyn.theta - ctrl.optTraj.theta;
        thetaDotErr = dyn.thetaDot - ctrl.optTraj.thetaDot;
        s = thetaErr + ctrl.lambda*thetaDotErr;
        subplot(224), plot(dyn.t,s), xlabel('Time [s]'), ylabel('Sliding Surface'), grid
end

switch params.runParams.drive
    case 'DHD'
        % Plot energy loss per event
        figure, yyaxis right
        plot(eval.switchTimes(1:end-1),eval.loss,'*'), ylabel('Valve Loss [J]')
        yyaxis left
        plot(dyn.t,dyn.u), xlabel('Time [s]'), ylabel('Control Input [Nm]'), grid

        % Plot cummulative valve loss over time
        figure, plot(eval.switchTimes(1:end-1),cumsum(eval.loss),...
            dyn.t,eval.mechEnergy), legend('Switching Loss','Absorbed Energy')
        xlabel('Time [s]'), ylabel('Cummulative Energy [J]'), grid
end

% Plot cylinder velocity (add force later)
theta = dyn.theta;
thetaDot = dyn.thetaDot;
dLdt = params.hyd.dLdt(theta,thetaDot);
figure, plot(dyn.t,dLdt), xlabel('Time [s]'), ylabel('Cylinder Velocity [m/s]')
figure, plot(dyn.t,theta*180/pi), xlabel('Time [s]'), ylabel('Flap Position [deg]')


% Compare average powers
% For the EHA, show the losses as if the convex map was used
eha_str_loss = '';
eha_str_pow = '';
eha_str_RGP = '';
switch params.runParams.drive
    case 'EHA'
        eha_str_loss = [' (Estimated ', num2str(eval.aveLossHat/1e3,3), ' kW)'];
        eval.aveElecPowHat = eval.aveMechPow - eval.aveLossHat;
        eha_str_pow = [' (Estimated ', num2str(eval.aveElecPowHat/1e3,3), ' kW)'];
        eval.elecRGP_hat = eval.aveElecPowHat/ctrl.optTraj.avePow;
        eha_str_RGP = [' (Estimated ', num2str(eval.elecRGP_hat,2), ')'];
end

% Display the average powers
disp([params.runParams.drive,' with ' ,params.runParams.controller])
disp(['Optimal Power = ', num2str(ctrl.optTraj.avePow/1e3,3), 'kW']);
disp(['Mechanical Power = ', num2str(eval.aveMechPow/1e3,3), ' kW']);
disp(['Drivetrain Loss = ', num2str(eval.aveLoss/1e3,3), ' kW ',eha_str_loss]);
disp(['Electric Power = ', num2str(eval.aveElecPow/1e3,3), ' kW ',eha_str_pow]);
disp('Relative Generated Power:')
eval.mechRGP = eval.aveMechPow/ctrl.optTraj.avePow;
eval.elecRGP = eval.aveElecPow/ctrl.optTraj.avePow;
    disp(['      Mechanical: ', num2str(eval.mechRGP,2)])
    disp(['      Electrical: ',num2str(eval.elecRGP,2),eha_str_RGP]);