# FIDIMAG:
from micro import Sim, FDMesh, UniformExchange, Demag, DMI
from pc import NEB_Sundials

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


# Simulation Function
def relax_neb(k, maxst, simname, init_im, interp, save_every=10000):
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
    sim.Ms = cylinder

    sim.add(UniformExchange(A=A))

    # Bulk DMI
    sim.add(DMI(D=D))

    # Demagnetization energy
    # sim.add(Demag())

    # Define many initial states close to one extreme. We want to check
    # if the images in the last step, are placed mostly in equally positions
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
              stopping_dmdt=1)

# Initial images and interpolation between them
init_im = [np.load('./relaxation/sk_up.npy'),
           np.load('./relaxation/sk_down.npy')
           ]
# 8 interpolations in between
interp = [16]

# Define different ks for multiple simulations
krange = ['1e10']

# Timing variables
f = open('timings.dat', 'w')
t = Timer()

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
