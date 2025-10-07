
rho = 1024;
Vol = 297;
g = 9.81;
r_cob = 4;
r_cog = 5;
Khs = rho*Vol*g*r_cob - m*g*r_cog;
I = 5.025e6;
Iinf = 1.734e7;
x_timeDomain = [1.847, 1.116, -0.934, 7.431];
x_frequencyDomain = [0.7045, 0.3652, 1.4028];

% Transfer function
s = tf('s');
Khat = x_frequencyDomain(1)*s*1e7/(s^2+2*x_frequencyDomain(2)*x_frequencyDomain(3)*s+x_frequencyDomain(3)^2);
G = 1/((I+Iinf)*s^2 + Khat*s + Khs);
figure, bode(G), hold on

% State space
A = [0,            -Khs/(I+Iinf),0,        -1/(I+Iinf);...
    1,               0, 0,        0;...
x_timeDomain(3), 0, 0, -x_timeDomain(1)*1e-7;...
x_timeDomain(4), 0, 1*1e-7, -x_timeDomain(2)*1e-7]
B = [1/(I+Iinf);0;0;0]; C = [0,1,0,0];
sys = ss(A,B,C,0);
bode(sys), hold off, grid
