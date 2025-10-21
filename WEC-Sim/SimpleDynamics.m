
rho = 1024;
Vol = 297;
g = 9.81;
r_cob = 4;
r_cog = 5;
m = 127000;
Khs = rho*Vol*g*r_cob - m*g*r_cog;
I = 5.025e6;
Iinf = 1.734e7;
x_timeDomain = [1.9754, 1.1345, 7.6921];
x_frequencyDomain = [7.045, 0.3652, 1.4028];

% Transfer function
s = tf('s');
Khat = x_frequencyDomain(1)*s*1e7/(s^2+2*x_frequencyDomain(2)*x_frequencyDomain(3)*s+x_frequencyDomain(3)^2);
G = 1/((I+Iinf)*s^2 + Khat*s + Khs);
figure, bode(G), hold on

% State space
A = [0,            -Khs/(I+Iinf),0,        -1/(I+Iinf);...
    1,               0, 0,        0;...
    0,               0, 0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0, 1, -x_timeDomain(2)]
B = [1/(I+Iinf);0;0;0]; C = [0,1,0,0];
sys = ss(A,B,C,0);
bode(sys), hold off, grid, legend('SS','TF')

pole(G)
eig(A)
figure, step(G), hold on, step(sys), hold off

Ar = [0,-x_timeDomain(1);1,-x_timeDomain(2)];
Br = [0;x_timeDomain(3)]*1e7; Cr = [0,0,1];
figure, bode(Khat), hold on, bode(ss(Ar,Br,Cr,0)), hold off