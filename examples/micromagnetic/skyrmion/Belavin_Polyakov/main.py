import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, Demag, DMI, UniaxialAnisotropy
from fidimag.micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader

mu0 = 4 * np.pi * 1e-7


def init_m(pos):

    x, y = pos[0] - 500, pos[1] - 500

    if x**2 + y**2 < 60**2:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m(init_m)
    #sim.set_m((0,0.1,1))
    #sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    dmi = DMI(D=1.3e-3)
    sim.add(dmi)

    anis = UniaxialAnisotropy(-3.25e4, axis=(0, 0, 1))
    sim.add(anis)

    zeeman = Zeeman((0, 0, 6.014576e4))
    sim.add(zeeman, save_field=True)

    sim.relax(dt=1e-13, stopping_dmdt=0.5, max_steps=5000,
              save_m_steps=None, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def init_m_BP(pos):

    x, y = pos[0] - 500, pos[1] - 500
    l = 100

    r2 = x*x+y*y

    nx = 2*l*x/(r2+l*l)
    ny = 2*l*y/(r2+l*l)
    nz = (r2-l*l)/(r2+l*l)

    return (nx,ny,nz)


def relax_system_only_exchange(mesh):

    sim = Sim(mesh, name='relax_exchange_only')

    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m(init_m_BP)

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    sim.relax(dt=1e-13, stopping_dmdt=0.5, max_steps=5000,
              save_m_steps=None, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


if __name__ == '__main__':

    mesh = CuboidMesh(
        nx=501, ny=501, nz=1, dx=2.0, dy=2.0, dz=2.0, unit_length=1e-9, periodicity=(True, True, False))

    #relax_system(mesh)
    relax_system_only_exchange(mesh)

    # apply_field1(mesh)
    # deal_plot()
