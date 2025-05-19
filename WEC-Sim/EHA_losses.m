%% EHA Case

%% Define Force and velocity from WEC-sim
force = logsout.getElement('force').Values.Data;
velocity = logsout.getElement('velocity').Values.Data;
% electricPower = logsout.getElement('electricPower').Values.Data;
electricPower = zeros(size(time));
t = logsout.getElement('velocity').Values.Time;
dt = t(2) - t(1);


%% Calculate Actual Flow and delta P through/across variable displacement motor
%Calculate Flow
Q_Act = velocity*params.rodArea;
n = length(Q_Act);

% Pressure Differential (Pressure actuator - Pressure Rail)
deltaP = force/params.rodArea;


%% Fluid Properties
mu=(32e-6)*870;
B = 1.7e9;
rho = 870;

%% Define Pump Constants
% Angular Velocity
Wrpm = 2000; %revolutions per minute
w = Wrpm.*(2*pi/60); % radians per second

%% I am scaling the pump to be X times larger than the 107 cc b/c That will make -1<fracDisp<1
Scale = max(abs(Q_Act))/w*2*pi*1e6/107;
Scale = params.scaleHECM;
%Scale = 69.7309;  % Regular wave case, no codesign
%Scale = 55;  % Regular wave case, yes codesign
%Scale = 63.9467; % for irregular wave case, no codesign
%Scale = 130; % for irregular wave case, yes codesign
%% Manufacturer 107cc/rev
% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
    D = 107; % cc/rev
    d = (D*100^-3)/(2*pi); % m^3/rad 
% Torques Loss Constants
    Cf =  53.7e-3;
    Ch = 53.6;
    Cv = 23.5e3;
% Flow Loss Constants
    Cs = 4.26e-9;
    Cst = 0*1e-5;

%% Find fraction of total displacement
fracDisp = NaN(n,1);
T_Ideal = NaN(n,1);
TLoss = NaN(n,1);
T_Act = NaN(n,1);
P_out = NaN(n,1);
P_in = NaN(n,1);
QLoss = NaN(n,1);
P_leak = NaN(n,1);
P_mech = NaN(n,1);
P_L_elect = NaN(n,1);
P_out_sum = 0;
P_in_sum = 0;



% Need to check what quadrant this is in, because the equations differ
% QLoss and TLoss scale with the pump. Pump 10X larger => Qloss 10X larger 
for i = 1:n
    % Calculate fracDisp assuming it is a posative value
    fracDisp(i) = (Q_Act(i) + sign(deltaP(i))*Scale*abs(d*Cs*(deltaP(i))/mu) + sign(deltaP(i))*Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale-sign(deltaP(i))*Scale*abs(d*w*deltaP(i)/B));
        if fracDisp(i) <= 0 % if the assumption that fracdisp is + is incorrect, recalculate fracDisp assuming fracDisp is - 
            fracDisp(i) = (Q_Act(i) + sign(deltaP(i))*Scale*abs(d*Cs*(deltaP(i))/mu) + sign(deltaP(i))*Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale+sign(deltaP(i))*Scale*abs(d*w*deltaP(i)/B));
        end
    
    % Check if fracDisp was calculated correctly
    % Calculate the error due to the fractional displacement calculation
    % If I use the fracDisp I calculated to solve for Q_Act (which was given)
    % Do I get the same values?
    QLoss(i) = Scale*abs(d*Cs*(deltaP(i))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(i))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5);
    Q_Act_calc(i) = w*d*fracDisp(i)*Scale - sign(deltaP(i))*QLoss(i);
    
    
    T_Ideal(i) = deltaP(i)*d*fracDisp(i)*Scale;
    TLoss(i) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(i))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
    
    
    Q_Ideal(i) = fracDisp(i)*d*w*Scale;
    T_Act(i) = T_Ideal(i) + sign(w)*TLoss(i); % |T_Act| needs be < |T_Ideal|
    
    % Power out with 90% effiency
    if T_Act(i) < 0
        P_out(i) = .9*w*T_Act(i);
        P_L_elect(i) = -.1*w*T_Act(i); % the negative sign makes sure the loss is posative 
    else
        P_out(i) = w*T_Act(i)/.9;
        P_L_elect(i) =  (1/.9-1)*w*T_Act(i); % This loss accounts for energy that needs to come FROM the generator to the system.
    end
    P_in(i) = (force(i)*velocity(i)); 
    
    
    % check if the difference in work in and work out  (i.e. integral of the power between in and out) is the losses in the hydraulic and electrical components.
    P_leak(i) = abs(deltaP(i) *QLoss(i));
    P_mech(i) = w*TLoss(i);
   
end

EHA_Generator = max([ max(abs(P_out)) ,  max(abs(T_Act*w))])/1e3; % kW
EHA_Pump_size = D*Scale;

Work_In = sum(P_in)*dt;
Work_Out = sum(P_out)*dt;
e = Work_Out/Work_In;


B = 18; % m (This is the width of the oswec)
t_start = 50;
Energy_first_chunk = -sum(P_in(1:find(t==t_start)))*dt;
Total_energy_in = -sum(P_in)*dt;
ave_power_in = (Total_energy_in-Energy_first_chunk)/(t(end)-t_start); % Ave power for last 100 s
CW = ave_power_in/waves.power;
CWR_in = CW/B;

Energy_first_chunk = -sum(P_out(1:find(t==t_start)))*dt;
Total_energy_out = -sum(P_out)*dt;
ave_power_out = (Total_energy_out-Energy_first_chunk)/(t(end)-t_start); % Ave power for last 100 s
CW = ave_power_out/waves.power;
CWR_out = CW/B;

netEfficiency = Total_energy_out/Total_energy_in;

%disp(['CWR in: ' num2str(CWR_in)])
%disp(['CWR out: ' num2str(CWR_out)])
%disp(['Efficiency: ' num2str(e)])
%disp(['EHA pump size: ' num2str(round(EHA_Pump_size)) ' cc'])
%disp(['EHA generator size: ' num2str(round(EHA_Generator)) ' kW'])


%%
plotEfficiencyMap(Scale), hold on
scatter(fracDisp,deltaP/1e6,'r','x'), hold off
%%
Work_in = -cumsum(P_in)*dt;
Work_out = -cumsum(P_out)*dt;

figure()
plot(t,Work_in, t,Work_out);
legend('Work in', 'Work out');
xlabel('Time(s)');
ylabel('Work (J)');
title(['Efficiency = ',num2str(sum(P_out)/sum(P_in)*100),'%'])
mean_power_out = -P_out_sum/length(t);
max_power_out = max(-P_out);
mean_power_in = -P_in_sum/length(t);
max_power_in = max(-P_in);


figure
plot(t,fracDisp)
xlabel('Time [s]')
ylabel('Fractional Displacement')
grid on


%% Plot drive cycle force and velocity
figure, yyaxis left, plot(t,force/1e6), ylabel('Force [MN]'), xlabel('Time [s]'), ylim([-ceil(max(abs(force/1e6))) ceil(max(abs(force/1e6)))]), grid on
yyaxis right, plot(t,velocity), ylabel('Velocity [m/s]'), ylim([-ceil(max(abs(velocity))) ceil(max(abs(velocity)))]), grid on


% figure
% yyaxis left
% plot(t,-P_in/1000,t,-P_out/1000)
% xlabel('Time (s)')
% ylabel('Power (kW)')
% %ylim([-2e5 12e5])
% yyaxis right
% plot(t,v)
% legend('Power In','Power Out','Velocity','location','northwest')
% xlim([150 200])
% ylim([-2 12])
% ylabel('Velocity (m/s)')
% grid on

figure
plot(t,-P_in/1000,t,-P_out/1000,t,electricPower/1e3)
xlabel('Time [s]')
ylabel('Power [kW]')
legend('Power In','Power Out','location','northwest')
grid on

% figure()
% plot(t,T_Act,t,T_Ideal)
% ylabel('Torque (Nm)'), xlabel('Time (s)')
% legend('Actual Torque','Ideal Torque')
% 
% figure()
% plot(t,Q_Act,t,Q_Ideal)
% ylabel('(m^3/s)'), xlabel('Time (s)')
% legend('Actual flow','Ideal flow')

%figure(), plot(t,T_Act,t,deltaP*d*Scale.*fracDisp + sign(w).*TLoss), legend('True','sum'), ylabel('Torque'), title('Torque')

%figure(), plot(t,Q_Act,t,w*d*Scale.*fracDisp - sign(deltaP).*QLoss), legend('true','sum'), ylabel('Flow'), title('Flow')

%figure, plot(t,cumsum(-P_in)*dt)

%% Instantaneous efficiency
Hydraulic_in = (-P_in>=0).*(-P_in);
Hydraulic_out = -(-P_in<=0).*(-P_in);
Mechanical_in = -(-P_out<=0).*(-P_out);
Mechanical_out = (-P_out>=0).*(-P_out);
inst_e = -(Hydraulic_out-Mechanical_out)./(Hydraulic_in-Mechanical_in);
figure, plot(t,inst_e), ylabel('Instantaneous efficiency'), xlabel('Time (s)')

ave_Hydraulic_in = sum((-P_in>=0).*(-P_in))*dt;
ave_Hydraulic_out = sum(-(-P_in<=0).*(-P_in))*dt;
ave_Mechanical_in = sum(-(-P_out<=0).*(-P_out))*dt;
ave_Mechanical_out = sum((-P_out>=0).*(-P_out))*dt;
ave_inst_e = -(ave_Hydraulic_out-ave_Mechanical_out)/(ave_Hydraulic_in-ave_Mechanical_in)

perc_pos_power_in = ave_Hydraulic_in/(ave_Hydraulic_in+ave_Hydraulic_out)

function plotEfficiencyMap(Scale)
%% Pump Motor Efficiency Mapping

%% Fluid Properties

mu=(32*10^-6)*870;
B = 1.7*10^9;
rho = 870;

%% Define Pump Constants

%%

%% UW Madison
% % Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
%     D = 107; % cc/rev
%     d = (D*100^-3)/(2*pi); % m^3/rad
% % Torques Loss Constants
%     Cf = 4.8*10^-3;
%     Ch = 0;
%     Cv = 153*10^3;
% % Flow Loss Constants
%     Cs = 1.04*10^-9;
%     Cst = 1.20*10^-5;
    
%% Manufacturer 107cc/rev
% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
    D = 107; % cc/rev
    d = (D*100^-3)/(2*pi); % m^3/rad
% Torques Loss Constants
    Cf = 53.7*10^-3;
    Ch = 53.6;
    Cv = 23.5*10^3;
% Flow Loss Constants
    Cs = 4.26*10^-9;
    Cst = 0*10^-5;

%% Define Pump Parameters to Study
% (fractional displacement/angular velocity/pressure differential)

% Angular Velocity
Wrpm = 3000;
w = Wrpm.*(2*pi/60); % radians per second

% Fractional Displacement Range
fracDisp = linspace(-1,1,50);

% Pressure Differential (Pressure actuator - Pressure Rail)
deltaP = -40e6:1e6:40e6;
% deltaP = [.5 1]*35e6;
%% Find Torque and Flow Losses

QLoss = zeros(length(fracDisp),length(deltaP));
Q_Ideal = zeros(length(fracDisp),length(deltaP));
Q_Act = zeros(length(fracDisp),length(deltaP));

TLoss = zeros(length(fracDisp),length(deltaP));
T_Ideal = zeros(length(fracDisp),length(deltaP));
T_Act = zeros(length(fracDisp),length(deltaP));
Q_Efficiency = zeros(length(fracDisp),length(deltaP));
T_Efficiency = zeros(length(fracDisp),length(deltaP));

Quad = zeros(length(fracDisp),length(deltaP));


for i = 1:length(fracDisp)
    for j = 1:length(deltaP)
         
%         Q_Ideal(i,j) = w(i)*d*fracDisp;
%         T_Ideal(i,j) = deltaP(j)*d*fracDisp;
%         QLoss(i,j) = abs(d*Cs*(deltaP(j))/mu) + abs(fracDisp*d*w(i)*(deltaP(j))/B) + abs(d^(2/3)*Cst*(2*(deltaP(j))/rho)^.5);
%         TLoss(i,j) = abs(d*Cv*mu*w(i)) + abs(d*(deltaP(j))*Cf) + abs(fracDisp*Ch*w(i)^2*rho*d^(5/3)/2);
        
        Q_Ideal(i,j) = fracDisp(i)*d*w*Scale;
        T_Ideal(i,j) = deltaP(j)*d*fracDisp(i)*Scale;
        QLoss(i,j) = Scale*abs(d*Cs*(deltaP(j))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(j))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(j))/rho)^.5);
        TLoss(i,j) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(j))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
%         Q_Act(i,j) = Q_Ideal-sign(w(i)*deltaP(j))*QLoss(i,j)
%         T_Act(i,j) = T_Ideal+sign(w(i)*deltaP(j))*TLoss(i,j)
        
        if deltaP(j)>=0 && w >=0
            Quad(i,j) = 1;
            Q_Act(i,j) = Q_Ideal(i,j) - QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) + TLoss(i,j);
        
        elseif deltaP(j)<0 && w >=0
            Quad(i,j) = 4;
            Q_Act(i,j) = Q_Ideal(i,j) + QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) + TLoss(i,j);
        
        elseif deltaP(j)>=0 && w <0
            Quad(i,j) = 2;
            Q_Act(i,j) = Q_Ideal(i,j) - QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) - TLoss(i,j);
            
        else
            Quad(i,j) = 3;
            Q_Act(i,j) = Q_Ideal(i,j) + QLoss(i,j);            
            T_Act(i,j) = T_Ideal(i,j) - TLoss(i,j);
        end
        
        if abs(Q_Act(i,j)) >= abs(Q_Ideal(i,j))
            Q_Efficiency(i,j) = Q_Ideal(i,j)/Q_Act(i,j);
        else
            Q_Efficiency(i,j) = Q_Act(i,j)/Q_Ideal(i,j);
        end
        
       if abs(T_Act(i,j)) >= abs(T_Ideal(i,j))
           T_Efficiency(i,j) = T_Ideal(i,j)/T_Act(i,j);
       else
           T_Efficiency(i,j) = T_Act(i,j)/T_Ideal(i,j);
       end
        
     end
end
% QEfficiency = Q_Act./Q_Ideal;
Q_Efficiency = Q_Efficiency';
T_Efficiency = T_Efficiency';
Total_Efficiency = (Q_Efficiency.*T_Efficiency);

%% Plot Contours

%levels = 0:.05:.95;



% % Volumetric
% figure
% %levels = [.995 .99 .985 .98 .975 .97];
% levels = 0:.1:1;
% contour(w,(deltaP),Q_Efficiency,levels,'ShowText','on')
% grid on,
% ylabel('Delta Pressure [MPa]')
% xlabel('Angular Velocity [rad/s]')
% %set(gcf, 'Position', get(0, 'Screensize'));
% 
% % Torque
% figure
% %levels = .9:.01:.99;
% contour(w,(deltaP),T_Efficiency,levels,'ShowText','on')
% grid on
% ylabel('Delta Pressure [MPa]')
% xlabel('Angular Velocity [rad/s]')
% %set(gcf, 'Position', get(0, 'Screensize'));

% Overall
levels = .1:.05:1.05;
levels = .25:.25:1;
contour(fracDisp,deltaP/1e6,Total_Efficiency,levels,'ShowText','on')
grid on
ylabel('$\Delta$ P [MPa]','interpreter','latex')
xlabel('Fractional Displacement','interpreter','latex')
% set(gcf, 'Position', get(0, 'Screensize'));
end