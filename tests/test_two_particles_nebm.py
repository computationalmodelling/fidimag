from __future__ import print_function
import pytest

# FIDIMAG:
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, UniaxialAnisotropy
from fidimag.common.neb_cartesian import NEB_Sundials
from fidimag.common.neb_spherical import NEB_Sundials as NEB_Sundials_spherical
import numpy as np

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
mesh = CuboidMesh(nx=3,
                  ny=1,
                  nz=1,
                  dx=3, dy=3, dz=3,
                  unit_length=1e-9
                  )


# Simulation Function
def relax_neb(k, maxst, simname, init_im, interp, save_every=10000,
              coordinates='Cartesian'):
    """
    Execute a simulation with the NEB function of the FIDIMAG code, for an
    elongated particle (long cylinder)

    The simulations are made for a specific spring constant 'k' (a float),
    number of images 'init_im', interpolations between images 'interp'
    (an array) and a maximum of 'maxst' steps.
    'simname' is the name of the simulation, to distinguish the
    output files.

    --> vtks and npys are saved in files starting with the 'simname' string

    """

    # Prepare simulation
    # We define the cylinder with the Magnetisation function
    sim = Sim(mesh)
    sim.Ms = two_part

    #sim.add(UniformExchange(A=A))

    # Uniaxial anisotropy along x-axis
    sim.add(UniaxialAnisotropy(Kx, axis=(1, 0, 0)))

    # Define many initial states close to one extreme. We want to check
    # if the images in the last step, are placed mostly in equally positions
    init_images = init_im

    # Number of images between each state specified before (here we need only
    # two, one for the states between the initial and intermediate state
    # and another one for the images between the intermediate and final
    # states). Thus, the number of interpolations must always be
    # equal to 'the number of initial states specified', minus one.
    interpolations = interp

    if coordinates == 'Cartesian':
        neb = NEB_Sundials(sim,
                           init_images,
                           interpolations=interpolations,
                           spring=k,
                           name=simname)
    elif coordinates == 'Spherical':
        neb = NEB_Sundials_spherical(sim,
                                     init_images,
                                     interpolations=interpolations,
                                     spring=k,
                                     name=simname)

    neb.relax(max_steps=maxst,
              save_vtk_steps=save_every,
              save_npy_steps=save_every,
              stopping_dmdt=1e-2)


def mid_m(pos):
    if pos[0] > 4:
        return (0.5, 0, 0.2)
    else:
        return (-0.5, 0, 0.2)


# this test runs for over a minute
@pytest.mark.slow
def test_energy_barrier_2particles():
    # Initial images: we set here a rotation interpolating
    init_im = [(-1, 0, 0), mid_m, (1, 0, 0)]
    interp = [6, 6]

    # Define different ks for multiple simulations
    krange = ['1e8']

    for k in krange:
        # print 'Computing for k = {}'.format(k)
        relax_neb(float(k), 2000,
                  'neb_2particles_k{}_10-10int'.format(k),
                  init_im,
                  interp,
                  save_every=5000)

        # Relax the same system using spherical coordinates
        relax_neb(float(k), 2000,
                  'neb_2particles_k{}_10-10int_spherical'.format(k),
                  init_im,
                  interp,
                  save_every=5000,
                  coordinates='Spherical'
                  )

    # Get the energies from the last state
    data = np.loadtxt('neb_2particles_k1e8_10-10int_energy.ndt')[-1][1:]
    data_spherical = np.loadtxt(
        'neb_2particles_k1e8_10-10int_spherical_energy.ndt')[-1][1:]

    ebarrier = np.abs(np.max(data) - np.min(data)) / (1.602e-19)
    print(ebarrier)

    ebarrier_spherical = np.abs(np.max(data_spherical) -
                                np.min(data_spherical)) / (1.602e-19)
    print(ebarrier_spherical)

    # Analitically, the energy when a single particle rotates is:
    #   K V cos^2 theta
    # where theta is the angle of the direction of one particle with respect
    # to the anisotropy axis. For this case, the MEP is the rotation
    # of a single particle and then followed by the rotation of
    # the other one (asynchronous). Thus, the barrier for the
    # parameters is: (27e-27) * 1e5 / (1.602e-19) ~ 0.01685 eV
    # since the volume is 3x3x3 nm^3 and theta is 1 for the maximum value
    # of a single rotated particle

    assert ebarrier < 0.017
    assert ebarrier > 0.005

    assert ebarrier_spherical < 0.017
    assert ebarrier_spherical > 0.005


if __name__ == '__main__':
    test_energy_barrier_2particles()
