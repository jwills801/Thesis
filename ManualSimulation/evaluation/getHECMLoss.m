function eval = getHECMLoss(eval,params,dyn)

u_hyd = NaN(size(dyn.u));
for tInd = 1:length(dyn.t)
    % Hydraulic torque
    u_hyd(tInd) = params.hyd.Force2Torque(dyn.theta(tInd))*params.hyd.ptoForceOptions(dyn.uInd(tInd));
end
u_elec = dyn.u-u_hyd;

% Powers
powHyd = u_hyd.*dyn.thetaDot;
powElec = u_elec.*dyn.thetaDot;
powTot = dyn.u.*dyn.thetaDot;

E_hyd = cumtrapz(dyn.t,powHyd);
E_elec = cumtrapz(dyn.t,powElec);
E_tot = cumtrapz(dyn.t,powTot);

% Loss from HECM
eff = 0.85;
Ploss = (1-eff)*(-powElec);
PlossAfterRamp = Ploss(dyn.t > params.simu.rampTime);

% Energy Loss from main motor
Eloss_mm = (1-eff)*(-E_hyd(end)) * (E_hyd(end)<0);

eval.TotalLoss = eval.TotalLoss+trapz(dyn.t,Ploss) + Eloss_mm;
eval.TotalLossAfterRamp = eval.TotalLossAfterRamp + sum(PlossAfterRamp)*params.simu.dt + Eloss_mm;
eval.aveLoss = eval.TotalLossAfterRamp / (params.simu.finalTime - params.simu.rampTime);

% Plots
figure, plot(dyn.t,powHyd,dyn.t,powElec,dyn.t,powTot)
legend('Hydraulic','Electric','Total'), ylabel('Power [W]')

figure, plot(dyn.t,E_hyd,dyn.t,E_elec,dyn.t,E_tot)
legend('Hydraulic','Electric','Total'), ylabel('Energy [J]')

a=1;
end
