import numpy as np
import fidimag
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
import fidimag.common.constant as C


def init_m(r):
    """
    Initial magnetisation profile; a dot to obtain a skyrmion
    """
    R_sk = 2
    rho = np.sqrt(r[0] ** 2 + r[1] ** 2)
    if rho < R_sk:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


dx, dy, dz = 0.25, 0.25, 1
Lx, Ly, Lz = 20, 20, 1
mesh = CuboidMesh(nx=int(Lx/dx), ny=int(Ly/dy), nz=int(Lz/dz),
                  dx=dx, dy=dy, dz=dz,
                  x0=-Lx/2, y0=-Ly/2, z0=-Lz/2,
                  unit_length=1e-9)

sim = Sim(mesh)

Ms = 1.1e6
A = 2e-12
D = 3.9e-3
Ku = 2.5e6
Bz = 1.

sim.set_Ms(Ms)
sim.add(fidimag.micro.UniformExchange(A))
sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))
sim.add(fidimag.micro.Zeeman((0, 0, Bz / C.mu_0)))

# sim.add(fidimag.micro.DMI(D, dmi_type='interfacial'))

# For a C_n material, there is a kind of instability when one of the DM
# constants is larger than approx 0.7 times the other DM constant
sim.add(fidimag.micro.DMI([D, 0.6 * D], dmi_type='C_n'))

sim.set_m(init_m)

sim.driver.do_precession = False
sim.driver.alpha = 0.9

sim.relax()

# ----------------------------------------------------------------------------
# Plot the 2D system

import matplotlib.pyplot as plt

m = sim.spin.reshape(-1, 3)
plt.figure(figsize=(8, 8))
plt.quiver(mesh.coordinates[:, 0], mesh.coordinates[:, 1],
           m[:, 0], m[:, 1],
           m[:, 2], cmap='RdYlBu'
           )
plt.show()
