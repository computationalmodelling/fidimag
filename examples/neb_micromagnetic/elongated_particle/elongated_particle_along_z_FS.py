# FIDIMAG:
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, UniaxialAnisotropy, Demag
from fidimag.common.nebm_FS import NEBM_FS
import numpy as np
import logging
import matplotlib.pyplot as plt

# Material Parameters ---------------------------------------------------------

A = 10e-12
Kx = 3e5
Ms = 3.98e5

# Mesh ------------------------------------------------------------------------

# Define an elongated cylinder along the y direction
def cylinder(r, centre, radius):
    if (r[0] - centre[0]) ** 2. + (r[1] - centre[1]) ** 2 <= radius ** 2.:
        return Ms
    else:
        return 0

# Finite differences mesh
mesh = CuboidMesh(nx=6, ny=6, nz=35,
                  dx=2, dy=2, dz=2,
                  unit_length=1e-9)

centre = (np.max(mesh.coordinates[:, 0]) * 0.5,
          np.max(mesh.coordinates[:, 1]) * 0.5)


# Prepare simulation ----------------------------------------------------------

def elongated_part_sim():
    sim = Sim(mesh)
    sim.Ms = lambda r: cylinder(r, centre, 8)
    sim.add(UniformExchange(A=A))
    sim.add(UniaxialAnisotropy(Kx, axis=(0, 0, 1)))  # Anisotropy along y
    sim.add(Demag())

    return sim


# -----------------------------------------------------------------------------

init_im = [(0, 0, -1), (0, 0.9, 0.1), (0, 0, 1)]
interp = [10, 10]

sim = elongated_part_sim()
nebm = NEBM_FS(sim, init_im, interpolations=interp, name='neb_cylinder_z-axis_FS',
                interpolation_method='rotation', spring_constant=1e5)

# dt = integrator.stepsize means after every integrator step, the images
# are rescaled. We can run more integrator steps if we decrease the
# stepsize, e.g. dt=1e-3 and integrator.stepsize=1e-4
nebm.integrator.maxSteps = 10
nebm.integrator.run_for(10,
                        # save_vtks_every=save_every,
                        # save_npys_every=save_every,
                        )

bandEnergies = np.loadtxt('neb_cylinder_z-axis_FS_energy.ndt')[:, 1:]
bandDistances = np.loadtxt('neb_cylinder_z-axis_FS_dYs.ndt')[:, 1:]

# Get the Matplotlib logger
mpl_logger = logging.getLogger('matplotlib')
# Set the logging level to 'WARNING' or higher
mpl_logger.setLevel(logging.WARNING)
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.WARNING)

f, ax  = plt.subplots()
print(bandDistances[1])
dist = [1] + list(np.cumsum(bandDistances[1]))
ax.plot(dist, bandEnergies[0], 'o--')

st = -1
dist = [0] + list(np.cumsum(bandDistances[st]))
ax.plot(dist, bandEnergies[st], 'o-')
plt.show()

