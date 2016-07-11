from __future__ import print_function
import pytest

import numpy as np

# FIDIMAG libraries
from fidimag.common import CuboidMesh as Mesh
from fidimag.atomistic import Sim, UniformExchange, Anisotropy
import fidimag.common.constant as const

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

    if x < 450:
        return (1, 0, 0)
    elif 450 <= x < 550:
        return (0, 1, 1)
    else:
        return (-1, 0, 0)


def relax_system(mesh):

    # Only relaxation
    sim = Sim(mesh, name='relax')

    # Simulation parameters
    sim.set_tols(rtol=1e-8, atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5/mu0
    sim.mu_s = 1e-27/mu0
    sim.do_precession = False

    # The initial state passed as a function
    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    # Energies
    exch = UniformExchange(J=2e-20)
    sim.add(exch)

    anis = Anisotropy(0.01*2e-20, axis=(1,0,0))
    sim.add(anis)

    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Start relaxation and save the state in m0.npy
    sim.relax(dt=1e-14, stopping_dmdt=1e3, max_steps=5000,
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
    sim.set_tols(rtol=1e-12, atol=1e-12)
    sim.gamma = 2.211e5/mu0
    sim.mu_s = 1e-27/mu0
    sim.alpha = 0.05

    sim.set_m(np.load('m0.npy'))

    # Energies
    exch = UniformExchange(J=2e-20)
    sim.add(exch)

    anis = Anisotropy(0.01*2e-20, axis=(1,0,0))
    sim.add(anis)
    # dmi = DMI(D=8e-4)
    # sim.add(dmi)

    # Set the current in the x direction, in A / m
    # beta is the parameter in the STT torque
    sim.jx = -1e12
    sim.beta = 0.1

    # The simulation will run for x ns and save
    # 'snaps' snapshots of the system in the process
    ts = np.linspace(0, time * 1e-9, snaps)


    for t in ts:
        print('time', t)
        sim.run_until(t)
        #sim.save_vtk()
    np.save('m1.npy', sim.spin)


def test_stt_dw_atomistic():
    mesh = Mesh(nx=1000, ny=1, nz=1,
                dx=1.0, dy=1.0, dz=1.0,
                unit_length=1e-9)

    # Relax the initial state. It will save the last state
    # to the m0.npy file
    relax_system(mesh)

    m0 = np.load('m0.npy')
    m0.shape=(-1,3)
    mx0 = np.average(m0[:,0])

    excite_system(mesh)

    m1 = np.load('m1.npy')
    m1.shape=(-1,3)
    mx1 = np.average(m1[:,0])

    mu_s = (1e-27/mu0)
    v = 1e-27
    p = 0.5
    jx = 1e12
    alpha = 0.05
    beta = 0.1
    time = 0.1

    u = const.g_e * const.mu_B / (2 * const.c_e) * v/ mu_s * p * jx
    distance = (1+alpha*beta)/(1+alpha*alpha)*u*time 
    dmx = distance/1000*2

    print(mx0, mx1, dmx, abs(mx1-mx0-dmx))
    assert(mx1-mx0>0)
    assert(abs(mx1-mx0-dmx)<5e-5)

if __name__ == '__main__':
    test_stt_dw_atomistic()
