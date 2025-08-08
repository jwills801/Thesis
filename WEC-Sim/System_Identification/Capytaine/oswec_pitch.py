import os
os.environ["OMP_NUM_THREADS"] = "1" 
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt
import xarray as xr

## Define floating body of the flap
bem_file = os.getcwd() + os.path.sep + 'flap.dat'
mesh = cpt.load_mesh(mesh = bem_file)
# dof = cpt.rigid_body_dofs(rotation_center=(0, 0, -8.9))
dof = cpt.rigid_body_dofs()
cog = (0,0,-3.9)
flap = cpt.FloatingBody(mesh=mesh,
                        dofs=dof,
                        center_of_mass=cog)
# flap.keep_immersed_part()
# flap.show()

anim = flap.animate(motion={"Pitch": 0.2}, loop_duration=1.0)
# anim.run()

# # define the base (I don't think this is necessary)
# mesh = cpt.load_mesh(mesh = bem_file[1])
# dof = cpt.rigid_body_dofs()
# base = cpt.FloatingBody(mesh=mesh,
#                         dofs=dof,
#                         center_of_mass=bem_cg[1])

## Hydrostatics problem
hydrostatics = flap.compute_hydrostatics(rho=1023)
K = hydrostatics["hydrostatic_stiffness"]
print(K)

## Radiation problem
omega_range = np.linspace(9, 11, 3)
print('omega Range =', omega_range)
problems = [
    cpt.RadiationProblem(body=flap, radiating_dof=dof, omega=omega, water_depth=10.9,rho=1023)
    for dof in flap.dofs
    for omega in omega_range
]
solver = cpt.BEMSolver()
results = solver.solve_all(problems)
radiationData = cpt.assemble_dataset(results)
print('RADITATION')
print(radiationData.keys())
B = radiationData['radiation_damping'].sel(omega=omega_range[1])
print(B)
for i in np.linspace(0,5,6,dtype=int):
    print(i+1)
    print(B[i,i].values)
    print('______')



## Diffraction problem
problems = [
    cpt.DiffractionProblem(body=flap, omega=omega, water_depth=10.9,rho=1023)
    for dof in flap.dofs
    for omega in omega_range
]
solver = cpt.BEMSolver()
results = solver.solve_all(problems)
diffractionData = cpt.assemble_dataset(results)
print('DIFFRACTION')
print(diffractionData.keys())
K = diffractionData['hydrostatic_stiffness']
print('HYDROSTATIC STIFFNESS')
print(K)


## Plot the added mass of each dofs as a function of the frequency
# plot Added Mass
plt.figure()
for dof in flap.dofs:
    plt.plot(
        omega_range,
        radiationData['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )
plt.xlabel('omega [rad/s]')
plt.ylabel('added mass')
plt.legend()
plt.tight_layout()
# plt.show()

# plot Radiation Damping
plt.figure()
for dof in flap.dofs:
    plt.plot(
        omega_range,
        radiationData['radiation_damping'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )
plt.xlabel('omega [rad/s]')
plt.ylabel('radiation damping')
plt.legend()
plt.tight_layout()
# plt.show()

# plot Excitation Force
fig, axs = plt.subplots(2)
plt.suptitle('Excitation Force')
for dof in flap.dofs:
    axs[0].plot(
        omega_range,
        diffractionData['excitation_force'].sel(influenced_dof=dof).real,
        label=dof,
        marker='o',
    )
    axs[1].plot(
        omega_range,
        diffractionData['excitation_force'].sel(influenced_dof=dof).imag,
        label=dof,
        marker='o',
    )
axs[1].set_xlabel('omega [rad/s]')
axs[0].set_ylabel('real')
axs[1].set_ylabel('imag')
axs[0].legend()
plt.show()


