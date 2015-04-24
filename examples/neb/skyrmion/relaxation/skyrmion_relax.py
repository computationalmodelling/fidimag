# FIDIMAG:
from micro import Sim, FDMesh, UniformExchange, Demag, DMI
import numpy as np


# Material Parameters for FeGe
A = 8.78e-12
D = 1.58e-3
Ms = 3.84e5

radius = 25


# MESH
def cylinder(pos):

    # Relative position
    x, y = pos[0] - radius, pos[1] - radius

    if x ** 2 + y ** 2 < radius ** 2:
        return Ms
    else:
        return 0

# We will generate a 50 nm wide and 5nm thick disk
# Finite differences mesh
mesh = FDMesh(nx=25, ny=25, nz=5,
              dx=2, dy=2, dz=1,
              unit_length=1e-9
              )


# Initial magnetisation
def init_m(pos):

    x, y = pos[0] - radius, pos[1] - radius

    if x ** 2 + y ** 2 < radius ** 2:
        return (0, 0, 1)
    else:
        return (0, 0, -1)

# Prepare simulation
# We define the cylinder with the Magnetisation function
sim = Sim(mesh, name='skyrmion')
sim.do_procession = False
sim.alpha = 0.5
sim.Ms = cylinder

sim.set_m(init_m)

sim.add(UniformExchange(A=A))

# DMI
sim.add(DMI(D=D))

# Relax the system
sim.relax(dt=1e-12, stopping_dmdt=0.0001, max_steps=5000,
          save_m_steps=None,
          save_vtk_steps=None
          )

np.save('sk_up.npy', sim.spin)
sim.save_vtk()
