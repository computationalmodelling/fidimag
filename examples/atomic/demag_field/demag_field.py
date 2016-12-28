import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from pc import Sim, Demag
from fidimag.common import CuboidMesh
import fidimag.common.constant as const
from pc.save_vtk import SaveVTK


def mu_s(pos):
    x, y, z = pos
    x0, y0, r = 60 * 0.5, 60 * 0.5, 30 * 0.5

    if (x - x0)**2 + (y - y0)**2 <= r**2:
        return const.mu_s_1
    else:
        return 0


def init_m(pos):
    x, y, z = pos

    return (0, 0, 1)


def random_m(pos):
    return np.random.random(3) - 0.5


def plot_f(mesh, field, mu_s_inv):
    fig = plt.figure(figsize=(6, 4))

    field.shape = (3, -1)

    xs = range(121)
    fz = []
    for i in xs:
        j = mesh.index(i, 60, 0)
        if mu_s_inv[j] > 0:
            fz.append(0)
        else:
            fz.append(field[2][j])

    plt.plot(xs, fz)
    plt.savefig('f.pdf')


def relax_system():

    mesh = CuboidMesh(nx=121, ny=121, nz=1, dx=0.5, dy=0.5, unit_length=1e-9)

    sim = Sim(mesh, name='relax_skx')
    sim.set_default_options(gamma=const.gamma)

    sim.driver.alpha = 1.0

    sim.mu_s = mu_s

    sim.set_m(init_m)

    demag = Demag()

    sim.add(demag)

    field = demag.compute_field()

    print field

    vtk = SaveVTK(mesh, name='demag')
    vtk.save_vtk(field)

    sim.save_vtk()

    plot_f(mesh, field, sim._mu_s_inv)

    # np.save('m0.npy',sim.spin)


if __name__ == '__main__':

    np.random.seed(3)

    relax_system()
