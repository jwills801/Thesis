function out = makeSwitchLossMap(params)

% This code will make a tabulated grid of points 
% which can be interpolated between to find the energy loss
tic
%% Parameters
capArea = params.hyd.capArea;
hoseVolume = 0.1^2*pi*5; % m^3: Radius of hose squared times pi times length of the hose
stroke = params.hyd.stroke;
PR = params.hyd.pressureRails;

% Fluid properties
beta = 1.8e9; % Pa: bulk modulus
alpha = 0;% .03; % Percent entrained air

% Precompute hydraulic effort vectors
p0 = 101325; % Pa, atmospheric
HydraulicEffort = hydraulic_effort(beta,alpha);

% Valve constants
    % set valve constant k: Q = k*sqrt(deltaP)
valveConstant = 1*capArea/sqrt(2e6); % m^3/s/Pa: Speed of cylinder times area divided by sqrt of allowed pressure drop
valveDampingRatio = 0.7;
valveNaturalFrequency = 25*2*pi; % 25 Hz, 25*2*pi rad/s - takes about 20 ms to open
valveTF = tf(valveNaturalFrequency^2,[1,2*valveNaturalFrequency*valveDampingRatio,valveNaturalFrequency^2]);

% Time parameters
finalTime = 0.1;
dt = 1e-5;
tspan = 0:dt:finalTime;

% Set the delay indices to be used
    % These are amounts of indecies to be delayed
    % The amounds of time span from 0 to 0.02 seconds
maxdelay = round(.02/dt);
delayvals = linspace(1,maxdelay,11);

% Set the volume and velocity*Area parameters
velA_n = 11;
velA_vals = linspace(-1.5*capArea,1.5*capArea,velA_n);

vol_n = 11;
vol_vals = linspace(hoseVolume,hoseVolume+stroke*capArea,vol_n);

% Initialize matrices
Eloss = NaN(length(PR),length(PR),velA_n,vol_n); delaychosen = NaN(length(PR),length(PR),velA_n,vol_n);
Eloss_delay = NaN(size(delayvals));

waitbarhandel = waitbar(0,'Calculating Switching Losses');
for m = 1:length(vol_vals)
    waitbar((m-1)/vol_n,waitbarhandel)
    for i = 1:length(velA_vals)
        for j = 1:length(PR)
            for k = 1:length(PR)
                if j == k % Pressures are the same, no switch occurs
                    % Loss is just the open valve loss
                    Eloss(j,k,i,m) = (abs(velA_vals(i)))^3/(valveConstant^2)*finalTime;

                    % Choose the first delay index
                    delaychosen(j,k,i,m) = 1;
                else % The pressures are different - a switch occurs
                    
                    % Loop over delay values
                    % Initialize loss vector for each delay
                    Eloss_delay(:) = NaN;
                    for ii = 4% 1:length(delayvals)
                        delay = round(delayvals(ii));

                        % Get valve trajectories
                        xon = step(valveTF,tspan);
                        xoff = (1-xon);
                        xon = [zeros(delay,1);xon(1:end-delay)];

                        % The is no oscillation in the valve. once it opens, it stays open
                        % Once it closes it closes
                        xoff((find(xoff<0,1)):end) = 0;
                        xon((find(xon>1,1)):end) = 1;

                        % Calculate volume
                        vol = vol_vals(m) - velA_vals(i)*cumsum(dt*ones(1,length(tspan)-1));

                        % Simulate pressure dynamics
                        % Initialize vectors
                        PA = NaN(size(tspan')); Q_C = NaN(length(tspan)-1,1); Q_O = NaN(length(tspan)-1,1);
                        PA(1) = PR(j);
                        for kk=1:length(tspan)-1
                            Q_C(kk) = valveConstant*xoff(kk)*sign(PR(j)-PA(kk)).*sqrt(abs(PR(j)-PA(kk))); % Oriface equation for closing valve
                            Q_O(kk) = valveConstant*xon(kk)*sign(PR(k)-PA(kk)).*sqrt(abs(PR(k)-PA(kk))); % Oriface equation for opening valve
                            [beta_mix, ~, ~] = oil_comp(PA(kk)+p0, alpha, beta); % Effective bulk modulus
                            %beta_mix = beta_oil;
                            % PA(kk+1) = PA(kk) + beta_mix/vol_vals(m)*(Q_C(kk)+Q_O(kk)-velA_vals(i))*dt; % Compressibility -> expression for dPdt
                            PA(kk+1) = PA(kk) + beta_mix/vol(kk)*(Q_C(kk)+Q_O(kk)-velA_vals(i))*dt; % Compressibility -> expression for dPdt
                        end
                        % calculate losses
                        % Loss_h = Q_C.*(PR(j)-PA(1:length(tspan)-1));
                        % Loss_l = Q_O.*(PR(k)-PA(1:length(tspan)-1));

                        phi_h = interp1(HydraulicEffort.Pvals,HydraulicEffort.phivals,PR(j)+p0);
                        PA_sat = PA(1:length(tspan)-1); PA_sat(PA_sat<min(HydraulicEffort.Pvals)) = min(HydraulicEffort.Pvals); % Hydraulic effort function only goes so low in pressure
                        phi_A = interp1(HydraulicEffort.Pvals,HydraulicEffort.phivals,PA_sat(1:length(tspan)-1)+p0);
                        phi_l = interp1(HydraulicEffort.Pvals,HydraulicEffort.phivals,PR(k)+p0);
                        Loss_h = Q_C.*(phi_h-phi_A);
                        Loss_l = Q_O.*(phi_l-phi_A);
                        Eloss_delay(ii) =  sum(Loss_h+Loss_l)*dt;
                        if 0%m == 1 && i == 5 && j == 2 && k == 3 && ii == 4
                            figure, title(['Delay value: ',num2str(delay)])
                            subplot(121), plot(tspan,PA)
                            subplot(122), plot(tspan,xon,tspan,xoff)

                            figure, plot(tspan(1:end-1),Loss_h,tspan(1:end-1),Loss_l), legend('Closing','Opening')
                            min(Eloss_delay)
                        end
                    end
                    [Eloss(j,k,i,m), delaychosen(j,k,i,m)] = min(Eloss_delay);
                end
            end
        end
    end
end

close(waitbarhandel)
disp(['Making switching losses took ' num2str(toc) ' seconds'])

out = struct();
out.Eloss = Eloss;
out.valveConstant = valveConstant;
out.velA_vals = velA_vals;
out.vol_vals = vol_vals;
out.PR = PR;
out.finalTime = finalTime;
out.hoseVolume = hoseVolume;
end


function out = hydraulic_effort(beta_oil,alpha)
p0 = 101325; % Pa, atmospheric
pg_vals = logspace(5,8,10);
phi_vals = NaN(size(pg_vals));
if alpha == 0
    phi_vals = [0.0010,0.0022,0.0047, 0.0101,0.0218,0.0469,0.1013,0.2189,0.4749,1.0386]*1e8;
else
    waitbarhandel = waitbar(0,'Precomputing Hydraulic Effort');
    for i = 1:length(pg_vals)
        waitbar((i-1)/length(pg_vals),waitbarhandel);

        phi_vals(i) = phi(pg_vals(i),p0,alpha,beta_oil);
        %phi_vals(i)/pg_vals(i)
    end
    close(waitbarhandel)
end

% output results
out.Pvals = pg_vals;
out.phivals = phi_vals;

% % Plot results
% figure, semilogx(pg_vals+p0,phi_vals./pg_vals)
% xlabel('Pressure [Pa]'), ylabel('\phi(p_g,p_0)/p_g'), grid

% figure, semilogx(pg_vals+p0,phi_vals01./pg_vals,pg_vals+p0,phi_vals03./pg_vals,pg_vals+p0,phi_vals05./pg_vals)
% xlabel('Pressure [Pa]'), ylabel('\phi(p_g,p_0)/p_g'), grid
% legend('\alpha = 1%','\alpha = 3%','\alpha = 5%',Location='northwest')

end


function phi = phi(pg,p0,alpha,beta_oil)
dp = pg/1e2;
P_vals = p0:dp:(pg+p0);
phi = 0;
for i = 1:length(P_vals)
    phi = phi + exp(g(pg+p0,P_vals(i),alpha,beta_oil))*dp;
end
end

function g = g(p2,p1,alpha,beta_oil)
dp = 100;
P_vals = p1:dp:p2;
g = 0;
for i = 1:length(P_vals)
    g = g + dp/oil_comp(P_vals(i), alpha, beta_oil);
end
end

function [beta_mix, W, rv] = oil_comp(P, alpha, beta_oil)
% Computes the bulk modulus, compressibility energy, and volume ratio of
% the oil based on the entrained air. Used in the cam synthesis, MATLAB
% analysis of a single piston, and simulink models. 

    %OIL_COMP compressibility of oil-air mixture
    %   beta_mix : spot bulk modulus of mixture at P; [Pa]
    %   W : volumetric energy density [J/m^3] 
    %   rv : volume ratio
    %   Pt : Tank pressure
    %   P : Abs fluid pressure
    %   alpha : entrained air fraction
    %   beta_oil : bulk modulus of pure oil
    Pt = 101325; % Pa, atmospheric
    gamma = 1.4; % Air constant
    Pg = P-Pt;   % Gauge pressure
    r = P/Pt;    % Pressure ratio
    e = exp(Pg/beta_oil);
    rv = (1./e+alpha*r.^(-1/gamma))/(1+alpha);      %Volume ratio

    %Note: beta_mix as calculated by Hans Barkei
    beta_mix = rv*(1+alpha)./(1./(e*beta_oil)+alpha./(gamma*Pt.*r.^(1/gamma+1)));
    
    El=beta_oil.*(e-(1+Pg./beta_oil))+Pt.*(e-1);           %Energy in liquid
    Eg=P./(gamma-1).*(1-r.^(1/gamma-1));          %Energy in gas
    E=(El./e+alpha*Eg.*r.^(-1/gamma))./(1./e+alpha*r.^(-1/gamma));    %Energy in mixture
    W = E-(1./rv-1).*Pt;  %Compressibility (volumetric) energy density

    % Code to plot bulk modulus as function of entrained air and pressure
    % beta_oil = 1.8e9;
    % alpha = [0 .01 .03 .05];
    % p0 = 101325;
    % P = linspace(.1,40e6,50);
    % 
    % Beta_mix = NaN(length(alpha),length(P));
    % for i = 1:length(alpha)
    %     [Beta_mix(i,:), ~, ~] = oil_comp(p0, P, alpha(i), beta_oil);
    % end
    % figure, plot(P/1e6,Beta_mix/1e9), legend('0%','1%','3%','5%'), grid, xlabel('Pressure [MPa]'), ylabel('Beta [GPa]')

end

