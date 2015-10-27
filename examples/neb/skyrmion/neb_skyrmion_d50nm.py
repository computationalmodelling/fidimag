"""

Example of the Cartesian NEB method applied to a 50 nm wide and 5 nm thick FeGe
nanodisk.

The system has two global minima, represented by two energetically equivalent
confined skyrmions with their core pointing perpendicular to the disk plane

The model for this system only take sinto account the exchange and the DMI
energies

"""

# FIDIMAG:
from fidimag.micro import Sim, UniformExchange, Demag, DMI
from fidimag.common import CuboidMesh
from fidimag.pc import NEB_Sundials

import numpy as np

# For timing purposes --------------------------------
# The timings will be saved to a file with the
# NEB relaxations timings
import time


class Timer:

    def __init__(self):
        self.t1 = 0
        self.t2 = 0

    def start(self):
        self.t1 = time.clock()

    def end(self):
        self.t2 = time.clock()

    def interval(self):
        return self.t2 - self.t1

# ----------------------------------------------------

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


# Simulation Function
def relax_neb(k, maxst, simname, init_im, interp, save_every=10000):
    """
    Execute a simulation with the NEB function of the FIDIMAG code, for a
    nano disk

    The simulations are made for a specific spring constant 'k' (a float),
    number of images 'init_im', interpolations between images 'interp'
    (an array) and a maximum of 'maxst' steps.
    'simname' is the name of the simulation, to distinguish the
    output files.

    --> vtks and npys are saved in folders starting with the 'simname'

    """

    # Prepare simulation
    # We define the small cylinder with the Magnetisation function
    sim = Sim(mesh)
    sim.Ms = cylinder

    # Energies

    # Exchange
    sim.add(UniformExchange(A=A))

    # Bulk DMI --> This produces a Bloch DW - like skyrmion
    sim.add(DMI(D=D))

    # No Demag, but this could have some effect
    # Demagnetization energy
    # sim.add(Demag())

    # Initial images (npy files or functions)
    init_images = init_im

    # Number of images between each state specified before (here we need only
    # two, one for the states between the initial and intermediate state
    # and another one for the images between the intermediate and final
    # states). Thus, the number of interpolations must always be
    # equal to 'the number of initial states specified', minus one.
    interpolations = interp

    # Initiate the NEB algorithm driver
    neb = NEB_Sundials(sim,
                       init_images,
                       interpolations=interpolations,
                       spring=k,
                       name=simname)

    # Start the relaxation
    neb.relax(max_steps=maxst,
              save_vtk_steps=save_every,
              save_npy_steps=save_every,
              stopping_dmdt=1)


# ----------------- NEB DEFINITIONS -------------------------------------------

# Initial images. We get them from the fidimag relaxations (so we are sure to
# get the ground states
init_im = [np.load('./relaxation/sk_up.npy'),
           np.load('./relaxation/sk_down.npy')
           ]

# 16 interpolations in between
interp = [16]

# Define different ks for multiple simulations (we only use k = 1e10
# for now, but we could compare with different magnitudes and study
# the influence of this parameter in the results). Normally, k = 1e10
# works fine, but when the systems are stiffer, we can use up to k = 1e12
# Larger values could make the NEB very slow and overdrive the energy band,
# generating repeated transitions along the energy landscape.
# A very small value will make the images to accumulate in the extremes
# or around critical points
krange = ['1e10']

# Timing variables
f = open('timings.dat', 'w')
t = Timer()

# Initiate the NEB realaxation and compute the timing of this simulation
for k in krange:
    print 'Computing for k = {}'.format(k)
    t.start()
    relax_neb(float(k), 2000,
              'neb_nanodisk_d50nm_k{}'.format(k),
              init_im,
              interp,
              save_every=50,
              )
    t.end()
    f.write('k{} {} \n'.format(k, t.interval()))
    f.flush()
f.close()
