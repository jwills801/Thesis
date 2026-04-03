hydro = struct();
% hydro = readCAPYTAINE(hydro,'oswec_full.nc');
hydro = readCAPYTAINE(hydro,'oswec_new2.nc');
hydro = radiationIRF(hydro,40,[],[],[],[]);
hydro = radiationIRFSS(hydro,[],[]);
hydro = excitationIRF(hydro,75,[],[],[],[]);
writeBEMIOH5(hydro)
plotBEMIO(hydro)

%%
filename = 'oswec_new2.nc';
% information vector (gives names of variables)
info = ncinfo(filename);

% Get frequency
w = ncread(filename,'omega');

% Diffraction and Froude Krylov forces add up to excitation force
DF = ncread(filename,'diffraction_force');
FK = ncread(filename,'Froude_Krylov_force');
EF_test = DF + FK;

% Excitation Force
K_exc = ncread(filename,'excitation_force');

% 
F = squeeze(K_exc(1,1,:,1)) + 1i*squeeze(K_exc(1,1,:,2));
F_test = squeeze(EF_test(1,1,:,1)) + 1i*squeeze(EF_test(1,1,:,2));
figure
subplot(211)
semilogx(w,20*log10(abs(F)),w,20*log10(abs(F_test)))
grid, xlim([.16, 5])

subplot(212)
semilogx(w,angle(F)*180/pi,w,angle(F_test)*180/pi)
grid, xlim([.16, 5])

return
%%
B = squeeze(hydro.B(5,5,:));
A = squeeze(hydro.A(5,5,:));
w = transpose(hydro.w);
K = B + w.*(A-hydro.Ainf(5,5));

figure, plot(hydro.w,K)


