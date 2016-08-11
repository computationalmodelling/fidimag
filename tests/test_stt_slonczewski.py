from __future__ import print_function

import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman, Anisotropy
from fidimag.common import CuboidMesh, DataReader


def test_dynamic():

    mesh = CuboidMesh(nx=1, ny=1, nz=1)

    sim = Sim(mesh, name='dyn_spin', driver='llg_stt_cpp')
    # sim.set_options(rtol=1e-10,atol=1e-14)
    sim.driver.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m((0.8,0,-1))

    Kx = Anisotropy(Ku=-0.05, axis=(0, 0, 1), name='Kz')
    sim.add(Kx)

    sim.p = (0,0,1)

    sim.a_J = 0.0052
    sim.alpha = 0.1

    ts = np.linspace(0, 1200, 401)
    for t in ts:
        sim.run_until(t)


    mz = sim.spin[2]
    alpha, K, u = 0.1, 0.05, 0.0052
    print(mz, u/(2*alpha*K))

    #########################################################
    # The system used in this test can be solved analytically, which gives that mz = u/(2*alpha*K),
    # where K represents the easy-plane anisotropy.
    ###
    assert abs(mz - u/(2*alpha*K))/mz< 5e-4
    



if __name__ == '__main__':

    test_dynamic()
