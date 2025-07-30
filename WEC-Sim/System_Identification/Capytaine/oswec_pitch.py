import os
os.environ["OMP_NUM_THREADS"] = "1" 
import numpy as np
import capytaine as cpt

sphere = cpt.mesh_sphere(radius=1.0, center=(0, 0, -2), name="my sphere")
sphere.vertices[:10]  # First ten vertices.
sphere.faces_centers[5]  # Center of the sixth face (Python arrays start at 0).
sphere.faces_normals[5]  # Normal vector of the sixth face.
# sphere.show()

body = cpt.FloatingBody(mesh=sphere,
                        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -2)),
                        center_of_mass=(0, 0, -2))

anim = body.animate(motion={"Pitch": 0.1}, loop_duration=1.0)
# anim.run()


##
# Define OSWEC parameters ----------------------------------------------------#
bem_file = (os.getcwd() + os.path.sep + 'flap.dat',
            os.getcwd() + os.path.sep + 'base_shift.stl')       # mesh files, base_cut.stl, base.dat nemoh, .gdf wamit
bem_cg = ((0,0,-3.90),
          (0,0,-10.90))                                         # centers of gravity
bem_name = ('oswec_flap',
            'oswec_base')                                       # body names

bodies = []
bodies.append(cpt.FloatingBody.from_file(bem_file[0]))
bodies[0].center_of_mass = bem_cg[0]
bodies[0].keep_immersed_part()
axis = cpt.Axis(vector=(0, 1, 0), point=(0, 0, -8.9))
bodies[0].add_rotation_dof(axis=axis, name='hingePitch')

anim = bodies[0].animate(motion={"hingePitch": 0.1}, loop_duration=1.0)
# anim.run()


## Try a different way
mesh = cpt.load_mesh(mesh = bem_file[0])
dof = cpt.rigid_body_dofs(rotation_center=(0, 0, -8.9))

body = cpt.FloatingBody(mesh=mesh,
                        dofs=dof,
                        center_of_mass=bem_cg[0])

anim = body.animate(motion={"Pitch": 0.2}, loop_duration=1.0)
anim.run()

