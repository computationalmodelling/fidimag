"""

NEB Simulation for a testing model made of Fe-like atoms
arranged in a 1D spin chain. This system is described by
50 spins with interfacial DMI, sepparated by a distance of
0.27 nm

There are two ground states: a helicoid and an uniform state,
wich are used to define an energy band interpolating the
magnetisation configuration between them. To get the
ferromagnetic state, we use an uniaxia anisotropy
perpendicular to the chain direction

This script requires the two relaxed states, that can be
generated running the scripts inside the 'relaxation'
folder:
                                   generates
    ./relaxation/1D_spin_chain.py     -->      helicoid
    ./relaxation/1D_spin_chain_fm.py  -->      ferromagnetic

Material Parameters:
    J = 12 meV
    D = 2 meV
    ku = 0.5 meV (along +z)
    mu = 2 mu_B (magnetisation)


The system is relaxed up to 400 steps and two spring constants:
    10^10 and 10^11
which can be adjusted as desired

"""

# FIDIMAG:
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import DMI
from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Zeeman
from fidimag.atomistic import Anisotropy
import fidimag.common.constant as const

# Import NEB library
from fidimag.common.neb_cartesian import NEB_Sundials

import numpy as np

# For timing purposes --------------------------------
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

# 1D chain of 50 spins with a lattice constant of 0.27 A
mesh = CuboidMesh(nx=50,
              dx=0.27,
              unit_length=1e-9,
              )


# Simulation Function
def relax_neb(k, maxst, simname, init_im, interp, save_every=10000):
    """
    Execute a simulation with the NEB function of the FIDIMAG code

    The simulations are made for a specific spring constant 'k' (a float),
    number of images 'init_im', interpolations between images 'interp'
    (an array) and a maximum of 'maxst' steps.
    'simname' is the name of the simulation, to distinguish the
    output files.

    --> vtks and npys are saved in files starting with the 'simname' string

    """

    # Prepare simulation
    sim = Sim(mesh, name=simname)
    sim.gamma = const.gamma

    # magnetisation in units of Bohr's magneton
    sim.mu_s = 2. * const.mu_B

    # Exchange constant in Joules: E = Sum J_{ij} S_i S_j
    J = 12. * const.meV
    exch = UniformExchange(J)
    sim.add(exch)

    # DMI constant in Joules: E = Sum D_{ij} S_i x S_j
    D = 2. * const.meV
    dmi = DMI(D, dmi_type='interfacial')
    sim.add(dmi)

    # Anisotropy along +z axis
    ku = Anisotropy(Ku=0.5 * const.meV,
                    axis=[0, 0, 1],
                    name='ku')
    sim.add(ku)

    # Initial images
    init_images = init_im

    # Number of images between each state specified before (here we need only
    # two, one for the states between the initial and intermediate state
    # and another one for the images between the intermediate and final
    # states). Thus, the number of interpolations must always be
    # equal to 'the number of initial states specified', minus one.
    interpolations = interp

    neb = NEB_Sundials(sim,
                       init_images,
                       interpolations=interpolations,
                       spring=k,
                       name=simname)

    neb.relax(max_steps=maxst,
              save_vtk_steps=save_every,
              save_npy_steps=save_every,
              stopping_dmdt=1e-2)

# Initial images and interpolation between them
init_im = [np.load('./relaxation/relax_spin_chain.npy'),
           np.load('./relaxation/relax_spin_chain_fm.npy')
           ]

# 16 interpolations in between
interp = [16]

# Define different ks for multiple simulations
krange = ['1e10', '1e11']

# Timing variables
f = open('timings.dat', 'w')
t = Timer()

for k in krange:
    print 'Computing for k = {}'.format(k)
    t.start()
    relax_neb(float(k),
              400,
              'neb_1D-chain-spins_helicoid-fm_atomic_k{}'.format(k),
              init_im,
              interp,
              save_every=10000,
              )
    t.end()
    f.write('k{} {} \n'.format(k, t.interval()))
    f.flush()
f.close()
