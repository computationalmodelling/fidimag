
import numpy as np
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import DMI
from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Zeeman
from fidimag.atomistic import Anisotropy
import fidimag.common.constant as const

# Number of spins and lattice spacing (see mesh)
nx = 50
dx = 0.27

sim_name = 'relax_spin_chain_fm'


def relax_system():

    # 1D chain of 50 spins with a lattice constant of 0.27 A
    mesh = CuboidMesh(nx=nx,
                  dx=dx,
                  unit_length=1e-9,
                  # pbc='1d'
                  )

    # Initiate the simulation. PBCs are specified in the mesh
    sim = Sim(mesh, name=sim_name)
    sim.driver.gamma = const.gamma

    # magnetisation in units of Bohr's magneton
    sim.mu_s = 2. * const.mu_B

    # sim.set_options(gamma=const.gamma, k_B=const.k_B)

    # Initial magnetisation profile
    sim.set_m((0, 0, 1))

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

    # Faster convergence
    sim.driver.alpha = 0.5
    sim.do_precession = False

    sim.relax(dt=1e-13, stopping_dmdt=0.05,
              max_steps=700,
              save_m_steps=1000, save_vtk_steps=1000)

    # Save the last relaxed state
    np.save(sim_name + '.npy', sim.spin)


relax_system()
