function out = PIcontrol(params)
% PI control
H = freqresp(params.sys,2*pi/params.period);
Kp = real(1/H(1)');
Ki = -2*pi/5*imag(1/H(1)');

end