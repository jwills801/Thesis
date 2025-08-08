hydro_old = struct();
hydro_old = readCAPYTAINE(hydro_old,'oswec_full.nc');
[w,w_ind] = min(abs(hydro_old(1).w-10));
disp('Exisiting Damping Coefficients at 10 rad/s')
hydro_old(1).B(1:6,1:6,w_ind)

disp('___________________________________________________')

hydro_new = struct();
hydro_new = readCAPYTAINE(hydro_new,'oswec_new.nc');
[w,w_ind] = min(abs(hydro_new(1).w-10));
disp('New Damping Coefficients at 10 rad/s')
hydro_new(1).B(1:6,1:6,w_ind)

disp('___________________________________________________')

hydro_new2 = struct();
hydro_new2 = readCAPYTAINE(hydro_new2,'oswec_new2.nc');
[w,w_ind] = min(abs(hydro_new2(1).w-10));
disp('Newer Damping Coefficients at 10 rad/s')
hydro_new2(1).B(1:6,1:6,w_ind)

return
hydro = radiationIRF(hydro,40,[],[],[],[]);
hydro = radiationIRFSS(hydro,[],[]);
hydro = excitationIRF(hydro,75,[],[],[],[]);
writeBEMIOH5(hydro)
plotBEMIO(hydro)


