import os
os.environ["OMP_NUM_THREADS"] = "1" 
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt
import xarray as xr

## Define floating body of the flap
flap_file = os.getcwd() + os.path.sep + 'flap.dat'
# dof = cpt.rigid_body_dofs(rotation_center=(0, 0, -8.9))
# flap = cpt.FloatingBody(mesh=cpt.load_mesh(mesh = flap_file),
#                         dofs=dof,
#                         center_of_mass=(0,0,-3.9))

flap = cpt.FloatingBody(mesh=cpt.load_mesh(mesh=flap_file), 
                        center_of_mass=(0,0,-3.9))
flap.rotation_center = (0, 0, -8.9)
flap.add_rotation_dof(name='Pitch')
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
omega_range = np.linspace(0.1, 5, 100)
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

# plot Excitation Force
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




###

## --- Isolated Plot for Pitch Excitation Force ---
fig, axs = plt.subplots(2, figsize=(8, 6), sharex=True)
plt.suptitle('Excitation Force (Pitch Only)')

# Select only the Pitch DOF
# Capytaine names these Surge, Sway, Heave, Roll, Pitch, Yaw
dof_to_plot = 'Pitch' 

# Extract the complex data for the specific DOF and wave direction
pitch_force = dataset['excitation_force'].sel(influenced_dof=dof_to_plot).squeeze()
pitch_mag = 20*np.log10(np.abs(pitch_force))
pitch_phase = np.rad2deg(np.angle(pitch_force))

# Plot Real Component
axs[0].plot(
    omega_range,
    pitch_mag,
    label=f'{dof_to_plot} (Real)',
    color='blue',
    linewidth=1.5
)
axs[0].set_ylabel('Magnitude [dB]')
axs[0].grid(True)
axs[0].legend()
axs[0].set_xscale('log')

# Plot Imaginary Component
axs[1].plot(
    omega_range,
    pitch_phase,
    label=f'{dof_to_plot} (Imag)',
    color='red',
    linewidth=1.5
)
axs[1].set_ylabel('Phase [deg]')
axs[1].set_xlabel('Omega [rad/s]')
axs[1].grid(True)
axs[1].legend()
axs[1].set_xscale('log')

plt.tight_layout()
plt.show()

nan_count = np.isnan(pitch_force.values).sum()
print(f"Number of NaN values in Pitch excitation: {nan_count} out of {len(pitch_force)}")
