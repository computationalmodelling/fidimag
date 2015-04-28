import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim
from micro import FDMesh
from micro import UniformExchange, Demag
from micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader

mu0 = 4 * np.pi * 1e-7


def spatial_Ms(pos):
    x = pos[0]
    y = pos[1]

    if (x - 100)**2 + (y - 100)**2 <= 100**2:
        return 8.6e5
    else:
        return 0


def init_m(pos):

    x, y = pos[0] - 100, pos[1] - 100

    if x**2 + y**2 < 10**2:
        return (0, 0, -1)
    elif x > 0 and y > 0:
        return (-1, 1, 0)
    elif x < 0 and y > 0:
        return (-1, -1, 0)
    elif x < 0 and y < 0:
        return (1, -1, 0)
    else:
        return (1, 1, 0)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = spatial_Ms
    sim.do_procession = False

    sim.set_m(init_m)
    # sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    mT = 795.7747154594767

    ONE_DEGREE_PER_NS = 17453292.52

    sim.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
              save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def excite_system(mesh):

    sim = Sim(mesh, name='dyn')

    sim.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.01
    sim.gamma = 2.211e5
    sim.Ms = spatial_Ms

    # sim.set_m(init_m)
    sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    mT = 795.7747154594767
    sigma = 0.08e-9

    def gaussian_fun(t):

        return np.exp(-0.5 * (t / sigma)**2)

    zeeman = TimeZeeman((80 * mT, 0, 0), time_fun=gaussian_fun, name='hx')
    #zeeman = Zeeman((100*mT,0,0), name='hx')
    sim.add(zeeman, save_field=True)

    ts = np.linspace(0, 1e-9, 501)

    for t in ts:
        print 'time', t
        print 'length:', sim.spin_length()[0:200]
        sim.run_until(t)
        sim.save_vtk()

if __name__ == '__main__':

    mesh = FDMesh(nx=80, ny=80, nz=2, dx=2.5, dy=2.5, dz=5.0, unit_length=1e-9)

    # relax_system(mesh)

    excite_system(mesh)

    # apply_field1(mesh)
    # deal_plot()
