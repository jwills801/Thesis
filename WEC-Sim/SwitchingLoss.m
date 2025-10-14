% This code will make a tabulated grid of points 
% which can be interpolated between to find the energy loss
tic
%% Parameters
hoseVolume = 0.02^2*pi*5; % m^3: Radius of hose squared times pi times length of the hose
rodArea = .05^2*pi; % m^2: Radius squared times pi
capArea = 1.5*rodArea; % m^2: Area ratio times rod Area
PR = [0 5 10 20 30]*1e6;
nPR = length(PR);

% set valve constant k: Q = k*sqrt(deltaP)
valveConstant = 1*capArea/sqrt(2e6); % m^3/s/Pa: Speed of cylinder times area divided by sqrt of allowed pressure drop
valveDampingRatio = 0.7;
valveNaturalFrequency = 25*2*pi; % 50 Hz, 50*2*pi rad/s - takes about 10 ms to open
beta = 1.8e9; % Pa: bulk modulus 

Tau = 1/valveDampingRatio/valveNaturalFrequency; % Time constant
finalTime = round(5*Tau,2);
dt = 1e-5;
tspan = 0:dt:finalTime;
gs = tf(valveNaturalFrequency^2,[1,2*valveNaturalFrequency*valveDampingRatio,valveNaturalFrequency^2]);

% Set the delay indices to be used
maxdelay = round(2.5*Tau/dt);
delayvals = linspace(1,maxdelay,10);

% Set the volume and velocity*Area parameters
velA_n = 11;
velA_vals = linspace(-1.5*capArea,1.5*capArea,velA_n);

vol_n = 11;
vol_vals = linspace(hoseVolume,hoseVolume+5*capArea,vol_n);

Eloss = NaN(length(PR),length(PR),velA_n,vol_n); delaychosen = NaN(length(PR),length(PR),velA_n,vol_n);
PA = NaN(size(tspan')); Q_C = NaN(length(tspan)-1,1); Q_O = NaN(length(tspan)-1,1);
Eloss_delay = NaN(size(delayvals));
waitbarhandel = waitbar(0,'Calculating Switching Losses');
for velA_ind = 1:velA_n
    waitbar(velA_ind/velA_n,waitbarhandel);
    for vol_ind = 1:velA_n
        for oldPR = 1:length(PR)
            for newPR = 1:length(PR)
                if oldPR == newPR
                    Eloss(oldPR,newPR,velA_ind,vol_ind) = (abs(velA_vals(velA_ind)))^3/(valveConstant^2)*finalTime;
                    delaychosen(oldPR,newPR,velA_ind,vol_ind) = NaN;
                else
                    Eloss_delay(:) = NaN;
                    for ii = 1:length(delayvals)
                        delay = round(delayvals(ii));


                        xon = step(gs,tspan);
                        xoff = (1-xon);
                        xon = [zeros(delay,1);xon(1:end-delay)];

                        % The is no oscillation in the vavle. once it opens, it stays open
                        % Once it closes it closes
                        xoff((find(xoff<0,1)):end) = 0;
                        xon((find(xon>1,1)):end) = 1;

                        PA(1) = PR(oldPR);
                        for kk=1:length(tspan)-1
                            Q_C(kk) = valveConstant*xoff(kk)*sign(PR(oldPR)-PA(kk)).*sqrt(abs(PR(oldPR)-PA(kk))); % Oriface equation for closing valve
                            Q_O(kk) = valveConstant*xon(kk)*sign(PR(newPR)-PA(kk)).*sqrt(abs(PR(newPR)-PA(kk))); % Oriface equation for opening valve
                            PA(kk+1) = PA(kk) + beta/vol_vals(vol_ind)*(Q_C(kk)+Q_O(kk)-velA_vals(velA_ind))*dt; % Compressibility -> expression for dPdt
                        end
                        Loss_h = Q_C.*(PR(oldPR)-PA(1:length(tspan)-1));
                        Loss_l = Q_O.*(PR(newPR)-PA(1:length(tspan)-1));
                        Eloss_delay(ii) =  sum(Loss_h+Loss_l)*dt;
                    end
                    [Eloss(oldPR,newPR,velA_ind,vol_ind), delaychosen(oldPR,newPR,velA_ind,vol_ind)] = min(Eloss_delay);
                end
                % if t_step > finalTime
                %     Eloss(oldPR,newPR,velA_ind,vol_ind) = Eloss(oldPR,newPR,velA_ind,vol_ind)...
                %                                 +(abs(velA_vec(velA_ind)))^3/k_/k_*(t_step-finalTime);
                % end
            end
        end
    end
end
close(waitbarhandel)

%%
params = struct; params.Eloss = Eloss;
params.velA_vals = velA_vals; params.vol_vals = vol_vals; params.PR = PR;
coeffs = NaN(4,nPR,nPR);
figure(3), counter = 0;
for oldPR_ind = 1:nPR
    for newPR_ind = 1:nPR
        counter = counter+1;
        % coeffs(:,oldPR_ind,newPR_ind) = LeastSquares(oldPR_ind,newPR_ind,params);
        % coeffs(:,oldPR_ind,newPR_ind) = [1,0,0,1];
        % a = coeffs(1,oldPR_ind,newPR_ind); b = coeffs(2,oldPR_ind,newPR_ind); c = coeffs(3,oldPR_ind,newPR_ind); d = coeffs(4,oldPR_ind,newPR_ind);
        loss = squeeze(Eloss(oldPR_ind,newPR_ind,:,:))/1e3;
        for velA_ind = 1:velA_n
            for vol_ind = 1:vol_n
                delP = (PR(oldPR_ind)-PR(newPR_ind))/1e6;

                a = -(delP/300).^2+.01;
                b = -.6/30*delP;
                c = (delP/65).^2;
                d = .01;

                Q = velA_vals(velA_ind)*1e3;
                V = vol_vals(vol_ind)*1e3;
                Loss_hat(velA_ind,vol_ind) = a*delP^2 + b*Q + c*V + d*Q^2;
                Loss_hat(velA_ind,vol_ind) = a*delP^2 + b*Q + c*V + d*Q^2;
            end
        end
        error = (Loss_hat - loss);
        J(counter) = sum(error(:).^2);
        figure(3), subplot(nPR,nPR,counter), 
        surf(velA_vals,vol_vals,Loss_hat'), hold(gca,'on'), surf(velA_vals,vol_vals,loss'), hold(gca,'off')
        title(['Delta P = ',num2str((PR(oldPR_ind)-PR(newPR_ind))/1e6),'MPa']), xlabel('Flow [m^3/s]'), ylabel('Volume [m^3]'), zlabel('Loss [kJ]')
        delP_chosen(counter) = (PR(oldPR_ind)-PR(newPR_ind))/1e6;
        a_chosen(counter) = a;
        b_chosen(counter) = b;
        c_chosen(counter) = c;
        d_chosen(counter) = d;
    end
end
% coeffs(:,oldPR_ind,newPR_ind)

figure, plot(J), ylabel('Cost'), xlabel('Counter'), grid
%%
x = -30:30;
figure, plot(delP_chosen,a_chosen,'*',x,-(x/300).^2+.01), ylabel('a'), xlabel('delta P [MPa]'), grid
figure, plot(delP_chosen,b_chosen,'*',x,-.6/30*x), ylabel('b'), xlabel('delta P [MPa]'), grid
figure, plot(delP_chosen,c_chosen,'*',x,(x/65).^2), ylabel('c'), xlabel('delta P [MPa]'), grid
figure, plot(delP_chosen,d_chosen,'*',x,d*ones(size(x))), ylabel('d'), xlabel('delta P [MPa]'), grid


%% Least squares minimization
function coeffs = LeastSquares(oldRail,newRail,params)

params.loss = squeeze(params.Eloss(oldRail,newRail,:,:))/1e3;
params.delP = (params.PR(oldRail)-params.PR(newRail))/1e6;

a_vals = -.01:.001:.01;
b_vals = -.6:.1:.6;
c_vals = -.01:.01:.2;
d_vals = 0:.01:.02;
[a_mat,b_mat,c_mat,d_mat] = ndgrid(a_vals,b_vals,c_vals,d_vals);
a_vec = a_mat(:); b_vec = b_mat(:); c_vec = c_mat(:); d_vec = d_mat(:);
cost = NaN(size(a_vec));
waitbarhandel2 = waitbar(0,'Search space for coefficients');
for ind = 1:length(a_vec)
    waitbar(ind/length(a_vec),waitbarhandel2);
    cost(ind) = fun(a_vec(ind),b_vec(ind),c_vec(ind),d_vec(ind),params);
end
close(waitbarhandel2)
[M,I] = min(cost);
coeffs = [a_vec(I),b_vec(I),c_vec(I),d_vec(I)];
end


function J = fun(a,b,c,d,params)
loss = params.loss;
delP=params.delP;
velA_vals = params.velA_vals;
vol_vals = params.vol_vals;

for velA_ind = 1:length(velA_vals)
    for vol_ind = 1:length(vol_vals)
        Q = velA_vals(velA_ind)*1e3;
        V = vol_vals(vol_ind)*1e3;
        Loss_hat(velA_ind,vol_ind) = a*delP^2 + b*Q + c*V + d*Q^2;
    end
end
error = (Loss_hat - loss);
J = sum(error(:).^2);
end