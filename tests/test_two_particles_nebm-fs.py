import pytest
import logging

# FIDIMAG:
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, UniaxialAnisotropy, Demag
from fidimag.common.nebm_FS import NEBM_FS
import numpy as np
import matplotlib.pyplot as plt

# Material Parameters
# Parameters
A = 1e-12
Kx = 1e5
# Strong anisotropy
Ms = 3.8e5


"""
We will define two particles using a 4 sites mesh, letting the
sites in the middle as Ms = 0

"""


def two_part(pos):

    x = pos[0]

    if x > 6 or x < 3:
        return Ms
    else:
        return 0

# Finite differences mesh
mesh = CuboidMesh(nx=3, ny=1, nz=1, dx=3, dy=3, dz=3, unit_length=1e-9)

# Simulation Function
def relax_string(maxst, simname, init_im, interp, save_every=10000):
    """
    """

    # Prepare simulation
    # We define the cylinder with the Magnetisation function
    sim = Sim(mesh)
    sim.Ms = two_part

    # sim.add(UniformExchange(A=A))

    # Uniaxial anisotropy along x-axis
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))

    # Define many initial states close to one extreme. We want to check
    # if the images in the last step, are placed mostly in equally positions
    # init_images = init_im

    # Number of images between each state specified before (here we need only
    # two, one for the states between the initial and intermediate state
    # and another one for the images between the intermediate and final
    # states). Thus, the number of interpolations must always be
    # equal to 'the number of initial states specified', minus one.
    interpolations = interp

    nebm = NEBM_FS(sim, init_im, interpolations=interpolations, name=simname,
                   interpolation_method='rotation', spring_constant=1e5)

    # dt = integrator.stepsize means after every integrator step, the images
    # are rescaled. We can run more integrator steps if we decrease the
    # stepsize, e.g. dt=1e-3 and integrator.stepsize=1e-4
    nebm.integrator.maxSteps = maxst
    nebm.integrator.run_for(maxst,
                            # save_vtks_every=save_every,
                            # save_npys_every=save_every,
                            )

    return nebm


def mid_m(pos):
    if pos[0] > 4:
        return (0.5, 0, 0.2)
    else:
        return (-0.5, 0, 0.2)


def test_energy_barrier_2particles_string():
    # Initial images: we set here a rotation interpolating
    init_im = [(-1, 0, 0), (0.0, 0.2, 1.0), (1, 0, 0)]
    interp = [6, 6]

    barriers = []

    # Define different ks for multiple simulations
    # krange = ['1e8']

    string = relax_string(30,
                          'nebmfs_2particles_k1e8_10-10int',
                          init_im,
                          interp,
                          save_every=5000,
                          )

    bandEnergies = np.loadtxt('nebmfs_2particles_k1e8_10-10int_energy.ndt')
    bandDistances = np.loadtxt('nebmfs_2particles_k1e8_10-10int_dYs.ndt')
    barriers.append((np.max(bandEnergies[-1][1:]) - bandEnergies[-1][1]) / 1.602e-19)

    print('Energy barrier is:', barriers[-1])
    # assert np.abs(barriers[-1] - 0.016019) < 1e-5

    print(barriers)

    # Get the Matplotlib logger
    mpl_logger = logging.getLogger('matplotlib')
    # Set the logging level to 'WARNING' or higher
    mpl_logger.setLevel(logging.WARNING)
    pil_logger = logging.getLogger('PIL')
    pil_logger.setLevel(logging.WARNING)

    f, ax  = plt.subplots()
    print(bandDistances[1])
    dist = [0] + list(np.cumsum(bandDistances[1]))
    ax.plot(dist, bandEnergies[1], 'o-')

    st = -1
    dist = [0] + list(np.cumsum(bandDistances[st]))
    ax.plot(dist, bandEnergies[st], 'o-')
    plt.show()


if __name__ == '__main__':
    test_energy_barrier_2particles_string()
