import numpy as np
from pc import *

import time

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def plot_mxyz(ts, mx, my, mz, me, name):
    fig = plt.figure()
    plt.plot(ts, mx, '^-', label='mx')
    plt.plot(ts, my, '.-', label='my')
    plt.plot(ts, mz, 'o-', label='mz')
    plt.plot(ts, me, '>-', label='me')
    plt.legend()
    fig.savefig(name)


def temperature_test(T):
    ni = Nickel()
    ni.alpha = 0.1
    ni.D = 1.35e-26
    ni.mu_s = 2.16e-23
    ni.J = 6.16e-21

    (nx, ny, nz) = (24, 24, 24)

    mesh = CuboidMesh(nx=nx, ny=ny, nz=nz)

    sim = Sim(mesh, T=T, driver='sllg')
    sim.set_options(dt=1e-15, gamma=ni.gamma, k_B=ni.k_B)
    sim.mu_s = ni.mu_s

    exch = UniformExchange(ni.J)
    sim.add(exch)

    anis = Anisotropy(ni.D)
    sim.add(anis)

    # zeeman=Zeeman(1e2,(0,0,1))
    # sim.add(zeeman)

    # demag=Demag(mu_s=ni.mu_s)
    # sim.add(demag)

    sim.set_m((0, 0.6, 0.99))

    # vs=VisualSpin(sim)
    # vs.init()

    ts = np.linspace(0, 1e-11, 101)
    me = []
    mx = []
    my = []
    mz = []
    for t in ts:
        sim.run_until(t)
        # vs.update()
        av = sim.compute_average()
        mx.append(av[0])
        my.append(av[1])
        mz.append(av[2])
        me.append(np.sqrt(np.sum(av * av)))
        # time.sleep(0.001)
        print av, me[-1]
        sim.save_vtk()
    name = 'nx%d_ny%d_nz%d_T%g.png' % (nx, ny, nz, T)
    plot_mxyz(ts, mx, my, mz, me, name)


if __name__ == '__main__':
    temperature_test(100)
