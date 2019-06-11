from fidimag.extensions.c_clib import PyExchangeEnergy
from fidimag.extensions.c_clib import PyMicroSim
from fidimag.common import CuboidMesh

import time
import numpy as np

A = np.ones(9, dtype=np.double)
print(A)
Exch = PyExchangeEnergy(A)
# A[2] = 4
# Exch.printA()

nx, ny, nz = 3, 3, 1
n = nx * ny * nz
dx, dy, dz = 1, 1, 1
unit_length = 1.

mesh = CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

spin = np.ones(3 * n, dtype=np.double)
spin[::3] = 0
spin[1::3] = 0
spin[3 * 4:3 * 4 + 3] = [0, 1, 0]
Ms = np.ones(n)
Ms_inv = 1 / Ms

field = np.zeros(3 * n, dtype=np.double)
energy = np.zeros(n, dtype=np.double)

Exch.setup(mesh.nx, mesh.ny, mesh.nz, mesh.dx, mesh.dy, mesh.dz,
           mesh.unit_length,
           spin, Ms, Ms_inv,
           mesh.coordinates, mesh.neighbours,
           energy, field
           )

Exch.compute_field(0)
print(field)
print(energy)

sim_C = PyMicroSim()
sim_C.add_interaction(Exch)

# time.sleep(5)
