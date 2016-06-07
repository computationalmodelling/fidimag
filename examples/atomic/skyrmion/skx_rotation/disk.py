import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman
from fidimag.common import BatchTasks, DataReader, CuboidMesh
from fidimag.common.helper import plot_m
import fidimag.common.constant as const


def init_m(pos):
    x, y, z = pos

    x1,y1, r = 50, 70, 5 

    m1 = (0.05, 0.01, -1)
    m2 = (0, 0, 1)

    if (x - x1)**2 + (y - y1)**2 < r**2:
        return m1
    else:
        return m2

def spatial_mu(pos):

    x, y, z = pos
    x0, y0, r = 70, 70, 70
    if (x-x0)**2 + (y-y0)**2 < r**2:
        return const.mu_s_1
    else:
        return 0



def random_m(pos):
    return np.random.random(3) - 0.5

def spatial_H(pos):
    x, y, z = pos
    x0, y0, r0 = 70, 70, 70
    r = np.sqrt((x-x0)**2 + (y-y0)**2)
    J = 50 * const.k_B
    h = 0.04+0.02*(r0-r)/r0

    return (0,0, h * J / const.mu_s_1)

def relax_system_stage1():

    mesh = CuboidMesh(nx=140 , ny=140, nz=1)

    sim = Sim(mesh, name='relax', driver='llg')
    #sim.set_options(dt=1e-14, gamma=const.gamma, k_B=const.k_B)
    sim.alpha = 0.5
    sim.do_precession = False
    sim.gamma = const.gamma
    sim.mu_s = spatial_mu

    sim.set_m(init_m)

    J = 50 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.27 * J
    dmi = DMI(D)
    sim.add(dmi)

    zeeman = Zeeman(spatial_H)
    sim.add(zeeman)

    sim.relax(dt=1e-14, stopping_dmdt=1e10, max_steps=1000,
              save_m_steps=100, save_vtk_steps=10)

    np.save('skx.npy', sim.spin)
    plot_m(mesh, 'skx.npy', comp='z')


def relax_system_stage2():

    mesh = CuboidMesh(nx=140 , ny=140, nz=1)

    sim = Sim(mesh, name='dyn', driver='llg')
    sim.alpha = 0.1
    sim.do_precession = True
    sim.gamma = const.gamma
    sim.mu_s = spatial_mu

    sim.set_m(np.load('skx.npy'))

    J = 50 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.27 * J
    dmi = DMI(D)
    sim.add(dmi)

    zeeman = Zeeman(spatial_H)
    sim.add(zeeman)

    ts = np.linspace(0, 2e-9, 201)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        sim.save_m()
        print(t)


if __name__ == '__main__':
    # np.random.seed(125)
    relax_system_stage1()
    relax_system_stage2()

