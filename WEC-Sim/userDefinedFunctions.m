mechanicalEnergy = logsout.getElement('Energy').Values.Data;
% EHA_losses;
switch PTO
    case 'Passive Valving'
        % note that pressure is not optimized for these wave conditions yet
        figure, subplot(211), plot(time,railPower,time,mechanicalPower), xlabel('Time [s]'), ylabel('Power [W]')
            legend('Rail','Mechanical'), grid, xlim([0 20])

        railCummulativeEnergy = cumtrapz(time,railPower);
        mechanicalCummulativeEnergy = cumtrapz(time,mechanicalPower);
        subplot(212), plot(time,railCummulativeEnergy,time,mechanicalCummulativeEnergy), xlabel('Time [s]'), ylabel('Cumulative Energy [J]')
            legend('Rail','Mechanical'), grid, xlim([0 20])
end
