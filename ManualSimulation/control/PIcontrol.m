function u = PIcontrol(states,params)

H = freqresp(params.phys.sys , 2*pi/params.simu.peakPeriod);

Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');

thetaDot = states(1);
theta = states(2);
u_cont = -1*(Kp*thetaDot + Ki*theta);

% Discretize
[~,uInd] = min(abs(u_cont-params.hyd.ptoTorqueOptions(:)));
u = params.hyd.ptoTorqueOptions(uInd);
end