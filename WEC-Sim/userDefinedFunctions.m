mechanicalEnergy = logsout.getElement('Energy').Values.Data;
% EHA_losses;
switch PTO
    case ['Passive Valving','Active Valving','HHEA']
        figure, subplot(211), plot(time,railPower,time,mechanicalPower), xlabel('Time [s]'), ylabel('Power [W]')
        legend('Rail','Mechanical'), grid, xlim([0 20])

        railCummulativeEnergy = cumtrapz(time,railPower);
        mechanicalCummulativeEnergy = cumtrapz(time,mechanicalPower);
        subplot(212), plot(time,railCummulativeEnergy,time,mechanicalCummulativeEnergy), xlabel('Time [s]'), ylabel('Cumulative Energy [J]')
        legend('Rail','Mechanical'), grid, xlim([0 20])

        outputEnergy = railCummulativeEnergy(end);
    case 'EHA'
        figure, subplot(211), plot(time,electricPower,time,mechanicalPower), xlabel('Time [s]'), ylabel('Power [W]')
        legend('Electrical','Mechanical'), grid, xlim([0 20])

        electricCummulativeEnergy = cumtrapz(time,electricPower);
        mechanicalCummulativeEnergy = cumtrapz(time,mechanicalPower);
        subplot(212), plot(time,electricCummulativeEnergy,time,mechanicalCummulativeEnergy), xlabel('Time [s]'), ylabel('Cumulative Energy [J]')
        legend('Rail','Mechanical'), grid, xlim([0 20])
        
        outputEnergy = electricCummulativeEnergy(end);
            
end
