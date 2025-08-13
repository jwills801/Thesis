import os
os.environ["OMP_NUM_THREADS"] = "1" 
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt
import xarray as xr

## Define floating body of the flap
flap_file = os.getcwd() + os.path.sep + 'flap.dat'
dof = cpt.rigid_body_dofs(rotation_center=(0, 0, -8.9))
# dof = cpt.rigid_body_dofs()
flap = cpt.FloatingBody(mesh=cpt.load_mesh(mesh = flap_file),
                        dofs=dof,
                        center_of_mass=(0,0,-3.9))
# flap.add_rotation_dof(name='Pitch')
flap.keep_immersed_part()

# # define the base
# base_file = os.getcwd() + os.path.sep + 'base_shift.stl'
# base = cpt.FloatingBody(mesh=cpt.load_mesh(mesh = base_file),
#                         dofs=cpt.rigid_body_dofs(),
#                         center_of_mass=(0,0,-10.9))
# bodies = flap + base
# bodies.show()

anim = flap.animate(motion={"Pitch": 0.2}, loop_duration=1.0)
# anim.run()
# print(flap.dofs)
print(list(flap.dofs))
print(['Surge','Sway','Heave','Roll','Pitch','Yaw'])
print(list(flap.dofs)[4])
print(list(['Pitch']))

## Hydrostatics
hydrostatics = flap.compute_hydrostatics(rho=1023.0)
np.savetxt('KH.dat', hydrostatics['hydrostatic_stiffness'])          # Write hydrostatic stiffness to KH.dat file
print(hydrostatics['hydrostatic_stiffness'])

# Write the other hydrostatics data to Hydrostatics.dat file
vo = hydrostatics['disp_volume']
cb = hydrostatics['center_of_buoyancy']
cg = (0,0,-3.9)
f = open('Hydrostatics.dat','w')
for j in [0,1,2]:
    line =  f'XF = {cb[j]:7.3f} - XG = {cg[j]:7.3f} \n'
    f.write(line)
line = f'Displacement = {vo:E}'
f.write(line)
f.close()

## Radiation and Diffraction
omega_range = np.linspace(0.04, 20, 500)
test_matrix = xr.Dataset({
        "omega": omega_range,
        "wave_direction": [0.],
        "radiating_dof": list(flap.dofs),
        "water_depth": [10.9],
        "rho": [1025.],
        })
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, flap)
# print(dataset.keys())

data = cpt.io.xarray.separate_complex_values(dataset)
data["radiating_dof"] = data["radiating_dof"].astype(str)
data["influenced_dof"] = data["influenced_dof"].astype(str)
data.to_netcdf('oswec_new2.nc')


# ## Plot
# # plot Added Mass
# plt.figure()
# for dof in flap.dofs:
#     plt.plot(
#         omega_range,
#         dataset['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
#         label=dof,
#         marker='o',
#     )
# plt.xlabel('omega [rad/s]')
# plt.ylabel('added mass')
# plt.legend()
# plt.tight_layout()
# # plt.show()

# # plot Radiation Damping
# plt.figure()
# for dof in flap.dofs:
#     plt.plot(
#         omega_range,
#         dataset['radiation_damping'].sel(radiating_dof=dof, influenced_dof=dof),
#         label=dof,
#         marker='o',
#     )
# plt.xlabel('omega [rad/s]')
# plt.ylabel('radiation damping')
# plt.legend()
# plt.tight_layout()
# # plt.show()

# # plot Excitation Force
# fig, axs = plt.subplots(2)
# plt.suptitle('Excitation Force')
# for dof in flap.dofs:
#     axs[0].plot(
#         omega_range,
#         dataset['excitation_force'].sel(influenced_dof=dof).real,
#         label=dof,
#         marker='o',
#     )
#     axs[1].plot(
#         omega_range,
#         dataset['excitation_force'].sel(influenced_dof=dof).imag,
#         label=dof,
#         marker='o',
#     )
# axs[1].set_xlabel('omega [rad/s]')
# axs[0].set_ylabel('real')
# axs[1].set_ylabel('imag')
# axs[0].legend()
# plt.show()