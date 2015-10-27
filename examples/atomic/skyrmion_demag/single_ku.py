import numpy as np
from pccp.pc import Sim
from fidimag.common import CuboidMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange
from pccp.pc import Demag
from pccp.pc import Anisotropy
from pccp.pc import Constant

from pccp.pc.batch_task import BatchTasks

const = Constant()


def mu_s(pos):
    x, y, z = pos
    x0, y0, r = 60 * 0.5, 60 * 0.5, 60.5 * 0.5

    if (x - x0)**2 + (y - y0)**2 <= r**2:
        return const.mu_s_1
    else:
        return 0


def init_m(pos):
    x, y, z = pos

    x0, y0, r = 60 * 0.5, 60 * 0.5, 25 * 0.5

    m1 = (0.05, 0.01, -1)
    m2 = (0, 0, 1)

    if (x - x0)**2 + (y - y0)**2 < r**2:
        return m1
    else:
        return m2


def random_m(pos):
    return np.random.random(3) - 0.5


def relax_system():

    mesh = CuboidMesh(nx=121, ny=121, dx=0.5, dy=0.5, unit_length=1e-9)

    sim = Sim(mesh, name='relax_skx')
    sim.set_options(gamma=const.gamma, k_B=const.k_B)

    sim.alpha = 1.0

    sim.mu_s = mu_s

    sim.set_m(init_m)
    # sim.set_m(random_m)

    J = 50.0 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.09 * J
    dmi = DMI(D)
    sim.add(dmi)

    K = 5e-3 * J
    anis = Anisotropy(K, direction=(0, 0, 1), name='Ku')
    sim.add(anis)

    sim.add(Demag())

    ONE_DEGREE_PER_NS = 17453292.52

    sim.relax(dt=1e-12, stopping_dmdt=0.1 * ONE_DEGREE_PER_NS,
              max_steps=5000, save_m_steps=100, save_vtk_steps=100)

    # np.save('m0.npy',sim.spin)


if __name__ == '__main__':

    np.random.seed(3)

    relax_system()
