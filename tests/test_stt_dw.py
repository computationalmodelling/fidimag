from __future__ import print_function
import pytest

# Since we are testing this script for different versions, we need to rebuild
# Fidimag every time we call the tests
# import os
# os.system('make clean && make build')

import numpy as np

# FIDIMAG libraries
from fidimag.micro import Sim
try:
    from fidimag.common import CuboidMesh as Mesh
except:
    from fidimag.micro import FDMesh as Mesh

# The energies (we can use DMI in a future simulation)
from fidimag.micro import UniformExchange
from fidimag.micro import UniaxialAnisotropy
# from micro import DMI


mu0 = 4 * np.pi * 1e-7


def load_mz_npy(npy_file):

    m0_z = np.load(npy_file)
    m_magnitude = m0_z[0] ** 2 + m0_z[1] ** 2 + m0_z[2] ** 2

    if np.abs(m_magnitude - 1) < 1e-6:
        m0_z = m0_z.reshape(-1, 3)[:, 2]
    else:
        m0_z = m0_z.reshape(3, -1)[2]

    return m0_z


# Initial State, a rough DW in a 1D chain
def init_m(pos):

    x = pos[0]

    if x < 400:
        return (1, 0, 0)
    elif 400 <= x < 500:
        return (0, 1, 1)
    else:
        return (-1, 0, 0)


def relax_system(mesh):

    # Only relaxation
    sim = Sim(mesh, name='relax')

    # Simulation parameters
    sim.driver.set_tols(rtol=1e-8, atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    # The initial state passed as a function
    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    # Energies
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Start relaxation and save the state in m0.npy
    sim.relax(dt=1e-14, stopping_dmdt=0.00001, max_steps=5000,
              save_m_steps=None, save_vtk_steps=None)

    np.save('m0.npy', sim.spin)


# This function excites the system with a
# current in the x-direction using Spin Transfer Torque
# formalism
def excite_system(mesh, time=5, snaps=501):

    # Specify the stt dynamics in the simulation
    sim = Sim(mesh, name='dyn', driver='llg_stt')

    # Set the simulation parameters
    sim.driver.set_tols(rtol=1e-12, atol=1e-14)
    sim.alpha = 0.05
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5

    # Load the initial state from the npy file saved
    # in the realxation
    sim.set_m(np.load('m0.npy'))

    # Add the energies
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Set the current in the x direction, in A / m
    # beta is the parameter in the STT torque
    sim.driver.jx = -1e12
    sim.driver.beta = 1

    # The simulation will run for x ns and save
    # 'snaps' snapshots of the system in the process
    ts = np.linspace(0, time * 1e-9, snaps)

    for t in ts:
        print('time', t)
        sim.driver.run_until(t)
        sim.save_vtk()
        sim.save_m()


# this test runs for about 10 seconds
@pytest.mark.slow
def test_stt_dw():
    # We will crate a mesh with 1000 elements of 2x2x2 nm
    # in the x direction, and 1 along y and z
    # (so we have a 1D system)
    mesh = Mesh(nx=1000, ny=1, nz=1,
                dx=2, dy=2, dz=2.0,
                unit_length=1e-9)

    # Relax the initial state. It will save the last state
    # to the m0.npy file
    relax_system(mesh)

    m0_z = load_mz_npy('m0.npy')

    x = np.arange(len(m0_z))
    index_max = np.argmax(np.abs(m0_z))

    assert x[index_max] == 225
    assert np.abs(m0_z[index_max] - 0.705740362679) < 1e-8

    # Excite the system for 1.5 ns
    excite_system(mesh, 1.5, 151)

    m0_z = load_mz_npy('dyn_npys/m_100.npy')
    x = np.arange(len(m0_z))
    # Check that the DW is at the 242th x-position in the 100th snapshot
    print(x[np.argmax(np.abs(m0_z))])
    assert x[np.argmax(np.abs(m0_z))] == 242

    # Check that the DW is at the 251th x-position in the 150th snapshot
    m0_z = load_mz_npy('dyn_npys/m_150.npy')
    assert x[np.argmax(np.abs(m0_z))] == 251

if __name__ == '__main__':
    test_stt_dw()
