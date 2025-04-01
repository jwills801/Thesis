% https://www.youtube.com/watch?v=HkVLvh4tfZ4&t=2959s

t = logsout.get('position').Values.Time;
x = logsout.get('position').Values.Data;
v = logsout.get('velocity').Values.Data;
F_PTO = logsout.get('force').Values.Data;

% F_ex = output.bodies(1).forceExcitation;
% a=5;
% F_ex5 = F_ex(:,1)*a.*cos(x(:,5)) - F_ex(:,3)*a.*sin(x(:,5)) + F_ex(:,5) ;

% figure, yyaxis left, plot(t,F_ex5/1e6), ylabel('Excitation Force [MN]')
% yyaxis right, plot(t,v(:,5)), ylabel('Velocity [rad/s]'), xlabel('Time [s]'), xlim([100 125])

% F_PTO = myoutput.signals.values(:,17);
% figure, plot(t,F_PTO(1:length(t))), xlabel('Time [s]'), ylabel('Linear PTO Force [N]')

%% Part b)
%Estimate the autocorrelation sequence using the sample autocorrelation
% $$\hat R_X [m] = \frac{1}{N} \sum^{N-1}_{n=0} X[n]X^*[n - m]$$
N_samp = length(v); dt = t(2)-t(1);
[R_v, vlags] = xcorr(v); R_v = R_v/N_samp;
[R_F, Flags] = xcorr(F_PTO); R_F = R_F/N_samp;
[R_vF, vFlags] = xcorr(v,F_PTO); R_vF = R_vF/N_samp;
figure, stem(vlags,R_v)
xlabel('Lag'), ylabel('Output (v) Autocorrelation'),grid
figure, stem(vFlags,R_F)
xlabel('Lag'), ylabel('Input (F) Autocorrelation'),grid
figure, stem(vFlags,R_vF)
xlabel('Lag'), ylabel('crosscorrelation'),grid


%% Part c)
%Using the estimated correlation sequence, estimate the PSD of $X[n]$ by computing the Fourier transform of $\hat R[m]$.

w = linspace(0,4,1000)*dt; 
Sv = freqz(R_v,1,w);
SF = freqz(R_F,1,w);
SvF = freqz(R_vF,1,w);
figure, plot(w/pi,20*log10(abs(Sv)))
xlabel('Normalized Frequency'), ylabel('v Power Spectrum [dB]'),grid
figure, plot(w/pi,20*log10(abs(SF)))
xlabel('Normalized Frequency'), ylabel('F Power Spectrum [dB]'),grid
figure, plot(w/pi,20*log10(abs(SvF)))
xlabel('Normalized Frequency'), ylabel('vF Power Spectrum [dB]'),grid

H_squared = Sv./SF;
H = SvF./SF;
%gamma_squared = abs(SvF).^2./Sv./SF; % Coherence - one means estimate is accurate - zero means estimate is inaccurate
figure, plot(w/dt,20*log10(abs(H_squared)))
xlabel('Frequency [rad/s]'), ylabel('H^2 Power Spectrum [dB]'),grid
figure, plot(w/dt,20*log10(sqrt(abs(H_squared))))
xlabel('Frequency [rad/s]'), ylabel('\sqrt{H^2} Power Spectrum [dB]','Interpreter','tex'),grid
figure, plot(w/dt,20*log10(abs(H)))
xlabel('Frequency [rad/s]'), ylabel('H Power Spectrum [dB]'),grid
figure, plot(w/dt,angle(H)*180/pi)
xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'),grid

%[Phat,w] = freqz(R_v,1);
%figure, plot(w/pi,20*log10(abs(Phat)))
%xlabel('Normalized Frequency'), ylabel('Power Spectrum [dB]')

%%
figure, subplot(211), semilogx(w/dt,10*log10(abs(H)))
xlabel('Frequency [rad/s]'), ylabel('H Power Spectrum [dB]'),grid
subplot(212), semilogx(w/dt,angle(H)*180/pi)
xlabel('Frequency [rad/s]'), ylabel('Phase [deg]'),grid

figure, subplot(211), semilogx(w/dt,20*log10(abs(1./H)))
xlabel('Frequency [rad/s]'), ylabel('Z Power Spectrum [dB]'),grid
subplot(212), semilogx(w/dt,angle(1./H)*180/pi)
xlabel('Frequency [rad/s]'), ylabel('Z Phase [deg]'),grid


%% Compare to optimal PI gain results
% Results for high frequency waves
%w_vals = [ 0.5500    0.6556    0.7611    0.8667    0.9722    1.0778    1.1833 1.2889    1.3944    1.5000];
%Mag_dB = [  147.8157  148.7084  152.0858  155.3360  157.1040  157.8162  157.8739  157.4890  156.7452  155.7937]; %dB
%Phase = [  -78.9405  -69.8674  -65.8491  -60.6948  -51.2515  -40.4163  -30.5154 -22.1088  -15.7891  -11.9373]; %Degrees


% Results for near resonance
%w_vals =[ 0.2000    0.3000    0.4000    0.5000    0.6000];
%Mag_dB =[146.6624  131.4759  139.3198  146.2295  149.7998];
%Phase = [86.2871   61.1235  -75.1945  -79.8062  -77.2729];

% both results
w_vals =[ 0.2000    0.3000    0.4000    0.5000    0.6000 0.5500    0.6556    0.7611    0.8667    0.9722    1.0778    1.1833 1.2889    1.3944    1.5000];
Mag_dB =[146.6624  131.4759  139.3198  146.2295  149.7998 147.8157  148.7084  152.0858  155.3360  157.1040  157.8162  157.8739  157.4890  156.7452  155.7937];
Phase = [86.2871   61.1235  -75.1945  -79.8062  -77.2729 -78.9405  -69.8674  -65.8491  -60.6948  -51.2515  -40.4163  -30.5154 -22.1088  -15.7891  -11.9373];


figure
subplot(211), semilogx(w/dt,20*log10(abs(conj(1./H))),w_vals,Mag_dB,'*'), ylabel('Magnitude [dB]'), grid, xlim([.1 10])
subplot(212), semilogx(w/dt,angle(conj(1./H))*180/pi,w_vals,Phase,'*'), xlabel('Frequency [rad/s]'), ylabel('Phase [Degrees]'), grid, xlim([.1 10])


%%
wn = .3265;
zeta = .012;
k = 1.5e-8;
G = tf([k 0],[1 2*zeta*wn wn^2]);
[MAG,PHASE,W] = bode(G);


figure
subplot(211), semilogx(w/dt,20*log10(abs(conj(1./H))),w_vals,Mag_dB,'*',W,20*log10(squeeze(1./MAG))), ylabel('Magnitude [dB]'), grid, xlim([.1 10])
subplot(212), semilogx(w/dt,angle(conj(1./H))*180/pi,w_vals,Phase,'*',W,squeeze(PHASE)), xlabel('Frequency [rad/s]'), ylabel('Phase [Degrees]'), grid, xlim([.1 10])