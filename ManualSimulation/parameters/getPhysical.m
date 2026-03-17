function phys = getPhysical(~)


phys.rho = 1024;
phys.Vol = 297;
phys.g = 9.81;
phys.r_cob = 4;
phys.r_cog = 5;
phys.waterDepth = 8;
phys.m = 127000;

% hydrostatic stiffness
phys.Khs = phys.rho * phys.Vol * phys.g * phys.r_cob - ...
    phys.m * phys.g * phys.r_cog;

% Inertia
phys.I = 5.025e6;

% Radiation
phys.Iinf = 1.734e7;
x_timeDomain = [1.9754, 1.1345, 7.6921];

% State space
I_total = phys.I+phys.Iinf;
A = [0,            -phys.Khs/I_total, 0, -1/I_total;...
    1,               0,               0,        0;...
    0,               0,               0, -x_timeDomain(1);...
x_timeDomain(3)*1e7, 0,               1, -x_timeDomain(2)];

B = [1/I_total;0;0;0]; 
C = eye(4);
phys.sys = ss(A,B,C,0);

end