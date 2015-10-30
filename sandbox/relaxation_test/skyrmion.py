"""

Relaxation to get a skyrmion using Fidimag for a toy model made of
Fe-like atoms arranged in a square lattice. This system is described by 21 x 21
spins with interfacial DMI and a lattice constant of 5 Angstrom

A very strong magnetic field B is applied perpendicular to the sample.

Magnetic parameters:
    J = 10 meV      Exchange
    D = 6 meV       DMI
    B = 25 T        Magnetic Field
    mu_s = 2 mu_B   Magnetic moment


Created by David Cortes on Fri 02 Oct 2015
University of Southampton
Contact to: d.i.cortes@soton.ac.uk

"""

# Numpy utilities
import numpy as np

# FIDIMAG Simulation imports:
from fidimag.atomistic import Sim
from fidimag.common.cuboid_mesh import CuboidMesh
from fidimag.atomistic import DMI
from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Zeeman
from fidimag.atomistic import Constant

# Import physical constants from fidimag
const = Constant()


# Define an initial state for the magnetisation as a function of space
# This function must return a tuple with (mx, my, mz)

# Initial skyrmion radius
in_radius = 2


# Initial state To get a skyrmion. It is only a small region
# at the center of the square with inverted spins, which will
# be the skyrmion core
def init_m(pos):

    x, y = pos[0] - 5.5, pos[1] - 5.5

    if x ** 2 + y ** 2 < in_radius ** 2:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


# MESH --------------------------------------------------------------------
# This is a 21x21 spins in a square lattice with a lattice constant of 5
# angstrom and PBCs
mesh = CuboidMesh(nx=10, ny=10,
                  dx=1, dy=1,
                  unit_length=1e-9,
                  # periodicity=(True, True, False)
                  )

# -----------------------------------------------------------------------------
# SIMULATION ------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Initiate a simulation object. PBCs are specified in the mesh
sim = Sim(mesh, name='relax_sk')
# Use default gamma value
sim.gamma = const.gamma

# Magnetisation in units of Bohr's magneton
sim.mu_s = 2 * const.mu_B

# We could change the parameters using this option
# sim.set_options(gamma=const.gamma)

# Initial magnetisation profile from the function
sim.set_m((0, -0.8, -0.8))

# Exchange constant in Joules: E = Sum J_{ij} S_i S_j
# J = (10.0 / 1.) * const.meV
# exch = UniformExchange(J)
# sim.add(exch)

# DMI constant in Joules: E = Sum D_{ij} S_i x S_j
# D = (6 / 1.) * const.meV
# dmi = DMI(D, dmi_type='interfacial')
# sim.add(dmi)

# Zeeman constant in Tesla:
sim.add(Zeeman((0, 0, 25.)))

# Tune the damping for faster convergence
sim.alpha = 0.5
# Remove precession
sim.do_precession = False

sim.run_until(1e-11)

# Relax the system
# The last state is saved automatically and we also save every 100 steps
# We can tune the LLG parameters and stopping criteria if necessary
# For the skyrmion we reduce the tolerances since
# the system produces relatively large errors using the default values
# (although the energy fluctuations are unnoticeable)
sim.set_tols(rtol=1e-10, atol=1e-12)
sim.relax(dt=1e-13,
          stopping_dmdt=0.01,
          max_steps=5000,
          save_m_steps=100, save_vtk_steps=100)
