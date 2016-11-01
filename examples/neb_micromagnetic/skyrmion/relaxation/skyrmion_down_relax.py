"""

Script to generate a skyrmion pointing DOWN in a 50 nm wide
and 5 nm thick FeGe disk

The skyrmion is stabilised due to the finite system (border
effects) and without Demag, external magnetic fields or
anisotropies

"""

# FIDIMAG:
from fidimag.micro import Sim, UniformExchange, Demag, DMI
from fidimag.common import CuboidMesh
import numpy as np


# Material Parameters for FeGe
A = 8.78e-12
D = 1.58e-3
Ms = 3.84e5

# Radius of the nanodisk
radius = 25

# MESH
# We pass this function to the Ms property
# of the simulation, so spins outside the desired
# radius will have Ms = 0
def cylinder(pos):

    # Relative position
    x, y = pos[0] - radius, pos[1] - radius

    if x ** 2 + y ** 2 < radius ** 2:
        return Ms
    else:
        return 0

# We will generate a 50 nm wide and 5nm thick disk The finite difference
# elements are 2nmx2nm cubes along the disk plane and they have a thickness of
# 1 nm
# Finite differences mesh
mesh = CuboidMesh(nx=25, ny=25, nz=5,
              dx=2, dy=2, dz=1,
              unit_length=1e-9
              )


# Initial magnetisation profile to get the skyrmion
# We create a small core pointing in the -z direction
def init_m(pos):

    x, y = pos[0] - radius, pos[1] - radius

    if x ** 2 + y ** 2 < radius ** 2:
        return (0, 0, -1)
    else:
        return (0, 0, 1)

# Prepare simulation
# We define the cylinder with the Magnetisation function
sim = Sim(mesh, name='skyrmion_down')
sim.Ms = cylinder

# To get a faster relaxation, we tune the LLG equation parameters
sim.do_precession = False
sim.alpha = 0.5

# Initial magnetisation:
sim.set_m(init_m)

# Energies:

# Exchange
sim.add(UniformExchange(A=A))

# Bulk DMI
sim.add(DMI(D=D))

# Relax the system
sim.relax(dt=1e-12, stopping_dmdt=0.0001, max_steps=5000,
          save_m_steps=None,
          save_vtk_steps=None
          )

# Save the final state and a vtk file
np.save('sk_down.npy', sim.spin)
sim.save_vtk()
