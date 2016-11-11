import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, Demag, DMI
from fidimag.micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader

mu0 = 4 * np.pi * 1e-7


def init_m(pos):

    x, y = pos[0] - 500, pos[1] - 500

    if x**2 + y**2 < 50**2:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.driver.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    dmi = DMI(D=1e-3)
    sim.add(dmi)

    zeeman = Zeeman((0, 0, 2e4))
    sim.add(zeeman, save_field=True)

    sim.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
              save_m_steps=None, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


if __name__ == '__main__':

    mesh = CuboidMesh(
        nx=1001, ny=1001, nz=1, dx=1, dy=1, dz=2.0, unit_length=1e-9, periodicity=(True, True, False))

    relax_system(mesh)

    # apply_field1(mesh)
    # deal_plot()
