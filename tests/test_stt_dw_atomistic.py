from __future__ import print_function
import pytest

import numpy as np

# FIDIMAG libraries
from fidimag.common import CuboidMesh as Mesh
from fidimag.atomistic import Sim, UniformExchange, Anisotropy
import fidimag.common.constant as const

mu0 = 4 * np.pi * 1e-7


# Initial State, a rough DW in a 1D chain
def init_m(pos):

    z = pos[2]

    if z < 450:
        return (0, 0, 1)
    elif 450 <= z < 550:
        return (1, 1, 0)
    else:
        return (0, 0, -1)


def relax_system(mesh):

    # Only relaxation
    sim = Sim(mesh, name='relax')

    # Simulation parameters
    sim.driver.set_tols(rtol=1e-8, atol=1e-10)
    sim.alpha = 0.5
    sim.driver.gamma = 2.211e5 / mu0
    sim.mu_s = 1e-27 / mu0
    sim.driver.do_precession = False

    # The initial state passed as a function
    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    # Energies
    exch = UniformExchange(J=2e-20)
    sim.add(exch)

    anis = Anisotropy(0.01*2e-20, axis=(0, 0, 1))
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Start relaxation and save the state in m0.npy
    sim.relax(dt=1e-14, stopping_dmdt=1e4, max_steps=5000,
              save_m_steps=None, save_vtk_steps=None)

    np.save('m0.npy', sim.spin)
    #sim.save_vtk()


# This function excites the system with a
# current in the x-direction using Spin Transfer Torque
# formalism
def excite_system(mesh, time=0.1, snaps=11):

    # Specify the stt dynamics in the simulation
    sim = Sim(mesh, name='dyn', driver='llg_stt')

    # Set the simulation parameters
    sim.driver.set_tols(rtol=1e-12, atol=1e-12)
    sim.driver.gamma = 2.211e5 / mu0
    sim.mu_s = 1e-27 / mu0
    sim.alpha = 0.05

    sim.set_m(np.load('m0.npy'))

    # Energies
    exch = UniformExchange(J=2e-20)
    sim.add(exch)

    anis = Anisotropy(0.01*2e-20, axis=(0, 0, 1))
    sim.add(anis)
    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Set the current in the x direction, in A / m
    # beta is the parameter in the STT torque
    sim.driver.jz = -1e12
    sim.driver.beta = 0.1

    # The simulation will run for x ns and save
    # 'snaps' snapshots of the system in the process
    ts = np.linspace(0, time * 1e-9, snaps)

    for t in ts:
        print('time', t)
        sim.driver.run_until(t)
        #sim.save_vtk()
    np.save('m1.npy', sim.spin)

    print(np.load('m1.npy')[:100])

def test_stt_dw_atomistic():
    mesh = Mesh(nx=1, ny=1, nz=1000,
                dx=1.0, dy=1.0, dz=1.0,
                unit_length=1e-9)

    # Relax the initial state. It will save the last state
    # to the m0.npy file
    relax_system(mesh)

    m0 = np.load('m0.npy')
    m0.shape = (-1, 3)
    mx0 = np.average(m0[:, 2])

    excite_system(mesh)

    m1 = np.load('m1.npy')
    m1.shape = (-1, 3)
    mx1 = np.average(m1[:, 2])

    mu_s = (1e-27 / mu0)
    v = 1e-27
    p = 0.5
    jx = 1e12
    alpha = 0.05
    beta = 0.1
    time = 0.1

    u = const.g_e * const.mu_B / (2 * const.c_e) * v / mu_s * p * jx
    distance = (1 + alpha * beta) / (1 + alpha * alpha) * u * time
    dmx = distance / 1000 * 2

    print(mx0, mx1, dmx, abs(mx1 - mx0 - dmx))
    assert(mx1 - mx0 > 0)
    assert(abs(mx1 - mx0 - dmx) < 5e-5)

if __name__ == '__main__':
    test_stt_dw_atomistic()
