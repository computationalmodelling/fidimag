import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from fidimag.atomistic import Sim, FDMesh, DMI, UniformExchange, Zeeman
from fidimag.common import Constant, BatchTasks, DataReader
from fidimag.common.helper import plot_m


const = Constant()


def init_m(pos):
    x, y, z = pos

    x0, y0, r = 28, 16, 4

    x1 = x % x0
    y1 = y % y0

    m1 = (0.05, 0.01, -1)
    m2 = (0, 0, 1)

    if (x1 - r)**2 + (y1 - r)**2 < r**2:
        return m1
    elif (x1 - x0 / 2. - r)**2 + (y1 - y0 / 2. - r)**2 < r**2:
        return m1
    else:
        return m2


def random_m(pos):
    return np.random.random(3) - 0.5


def excite_system(T=0.1, H=0.15):

    mesh = FDMesh(nx=28 * 3, ny=16 * 5, nz=1, pbc='2d')

    sim = Sim(mesh, name='dyn', driver='sllg')
    sim.set_options(dt=1e-14, gamma=const.gamma, k_B=const.k_B)
    sim.alpha = 0.1
    sim.mu_s = const.mu_s_1

    sim.set_m(random_m)

    J = 50 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.5 * J
    dmi = DMI(D)
    sim.add(dmi)

    Hz = H * J / const.mu_s_1
    zeeman = Zeeman([0, 0, Hz])
    sim.add(zeeman)

    sim.T = J / const.k_B * T

    ts = np.linspace(0, 5e-11, 51)
    for t in ts:
        sim.run_until(t)
        # sim.save_vtk()

    np.save('m.npy', sim.spin)
    plot_m(mesh, 'm.npy', comp='z')


def compute_skx_num(T=0.1, H=0.15):

    data = DataReader('dyn.txt')

    return data['skx_num'][-1]


def plot_res(Ts, Hs, res):
    mesh_x, mesh_y = np.meshgrid(Ts, Hs)
    res.shape = (len(Ts), len(Hs))

    fig = plt.figure(figsize=(6, 4))

    im = plt.imshow(np.transpose(res), cmap=cm.RdBu, extent=[
                    Ts[0], Ts[-1], Hs[0], Hs[-1]], origin='lower', aspect='auto')
    plt.colorbar(im)

    plt.text(0.25, 0.4, 'FM')
    plt.text(0.05, 0.15, 'SkX')
    plt.text(0.1, 0.005, 'HL')

    plt.xlabel(r'$k_B T/J$')
    plt.ylabel('H')
    plt.title('Skyrmion number')
    plt.tight_layout()

    fig.savefig('res.pdf')


if __name__ == '__main__':
    # np.random.seed(125)

    Ts = np.linspace(0, 0.5, 11)
    Hs = np.linspace(0, 0.5, 11)
    print Ts, Hs

    task = BatchTasks(excite_system, 1, taskname='task')
    task.add_parameters('T', Ts)
    task.add_parameters('H', Hs)

    task.start()

    task.post_process(compute_skx_num)

    pars, res = task.get_res()

    print pars
    print res

    np.save('res.npy', res)

    plot_res(Ts, Hs, res)
