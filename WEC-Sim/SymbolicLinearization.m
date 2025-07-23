syms I_11 I_22 I_33 I_44 I_55 I_66 I_15 I_51 theta theta_dot theta_ddot r

x_ddot =  r*theta_ddot*cos(theta) - r*theta_dot^2*sin(theta);
z_ddot = -r*theta_ddot*sin(theta) - r*theta_dot^2*cos(theta);
I = [I_11  0    0    0   I_15  0;
      0   I_22  0    0    0    0;
      0    0   I_33  0    0    0;
      0    0    0   I_44  0    0;
     I_51  0    0    0   I_55  0;
      0    0    0    0    0   I_66];

forceInerita = I * [x_ddot; 0; z_ddot; 0; theta_ddot; 0];
torqueInerita = forceInerita(5) + forceInerita(1)*r*cos(theta) - forceInerita(3)*r*sin(theta)

torqueIneritaLinearized = subs(torqueInerita,[theta, theta_dot, theta_ddot],[0 0 0]) + ...
                          subs(diff(torqueInerita,theta),[theta, theta_dot, theta_ddot],[0 0 0])*(theta-0) + ...
                          subs(diff(torqueInerita,theta_dot),[theta, theta_dot, theta_ddot],[0 0 0])*(theta_dot-0) + ...
                          subs(diff(torqueInerita,theta_ddot),[theta, theta_dot, theta_ddot],[0 0 0])*(theta_ddot-0)

%% Plot difference between linear and nonlinear inertial torques
Time = 0:.01:100;
Theta = 20*pi/180 * sin(Time/20*2*pi);
Theta_dot = [0 diff(Theta)]/Time(2);
Theta_ddot = [0,0 diff(Theta_dot(2:end))]/Time(2);
figure, plot(Time,Theta,Time,Theta_dot,Time,Theta_ddot), grid, legend('$\theta$','$\dot \theta$','$\ddot \theta$','Interpreter','latex','FontSize', 14)
xlabel('Time [s]'), ylabel('[rad], [rad/s], or [rad/s/s]')

I_numerical = double(subs(I,[ I_11       I_33    I_55     I_15      I_51   I_22 I_44 I_66],...
                            [1272734 2.1067e+05 7752007 -1646946 -1648846    0    0    0]));
r = 5;
X_ddot = r*Theta_ddot.*cos(Theta) - r*Theta_dot.^2.*sin(Theta);
Z_ddot = -r*Theta_ddot.*sin(Theta) - r*Theta_dot.^2.*cos(Theta);
forceInerita_numerical = I_numerical * [X_ddot; zeros(1,length(Time)); Z_ddot; zeros(1,length(Time)); Theta_ddot; zeros(1,length(Time))];
torqueInerita_numerical = forceInerita_numerical(5,:) + forceInerita_numerical(1,:)*r.*cos(Theta) - forceInerita_numerical(3,:)*r.*sin(Theta);

torqueIneritaLinearized_numerical = Theta_ddot*(7752007 + -1648846*r + r*(-1646946 + 1272734*r));

figure, plot(Time,torqueInerita_numerical,Time,torqueIneritaLinearized_numerical), grid, xlabel('Time [s]'), ylabel('Inertial Torque [Nm]')
legend('Non-linear','Linearized','FontSize', 11), ylim([-1e6 1e6])

%% Buoyancy and Restoring Force
K = zeros(6,6); K(3,3) = 317844; K(5,5) = -1933551;
m = 127000; g= 9.81; rho = 1000; V = 297.3760;
forceRestoring_numerical = zeros(6,length(Time));
    forceRestoring_numerical(3,:) = K(3,3)*(r*cos(Theta)-r) + g*(m-rho*V);
    forceRestoring_numerical(5,:) = K(5,5)*Theta;
torqueRestoring_numerical = forceRestoring_numerical(5,:) + forceRestoring_numerical(1,:)*r.*cos(Theta) - forceRestoring_numerical(3,:)*r.*sin(Theta);
torqueRestoringLinearized_numerical = (K(5,5) - r*g*(m-rho*V))*Theta;
figure, plot(Time,torqueRestoring_numerical,Time,torqueRestoringLinearized_numerical), grid, xlabel('Time [s]'), ylabel('Restoring Torque [Nm]')
legend('Non-linear','Linearized','FontSize', 11)%, ylim([-1e6 1e6])
