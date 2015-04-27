import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim
from micro import FDMesh
from micro import UniformExchange, Demag
from micro import Zeeman, TimeZeeman
from pc import DataReader

mu0 = 4 * np.pi * 1e-7


def init_m(pos):

    x = pos[0]

    if x <= 2:
        return (1, 0, 0)
    elif x >= 4:
        return (0, 0, 1)
    else:
        return (0, 1, 0)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-10, atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.0e5
    sim.do_procession = False

    sim.set_m((1, 0.25, 0.1))
    # sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    sim.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
              save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def apply_field1(mesh):

    sim = Sim(mesh, name='dyn')

    sim.set_tols(rtol=1e-10, atol=1e-10)
    sim.alpha = 0.02
    sim.gamma = 2.211e5
    sim.Ms = 8.0e5

    sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    mT = 0.001 / mu0
    print mT

    zeeman = Zeeman([-24.6 * mT, 4.3 * mT, 0], name='H')
    sim.add(zeeman, save_field=True)

    ts = np.linspace(0, 1e-9, 201)
    for t in ts:
        sim.run_until(t)
        print 'sim t=%g' % t


def deal_plot():
    data = DataReader('dyn.txt')
    ts = data['time'] * 1e9
    mx = data['m_x']
    my = data['m_y']
    mz = data['m_z']

    data2 = np.loadtxt('oommf.txt')

    ts2 = data2[:, 0] * 1e9
    mx2 = data2[:, 1]
    my2 = data2[:, 2]
    mz2 = data2[:, 3]

    plt.plot(ts, mx, '--', label='m_fidimag', dashes=(2, 2))
    plt.plot(ts, my, '--', label='', dashes=(2, 2))
    plt.plot(ts, mz, '--', label='', dashes=(2, 2))
    plt.plot(ts2, mx2, '--', label='m_oommf')
    plt.plot(ts2, my2, '--', label='')
    plt.plot(ts2, mz2, '--', label='')

    plt.title('std4')
    plt.legend()
    #plt.xlim([0, 0.012])
    #plt.ylim([-5, 100])
    plt.xlabel(r'Ts (ns)')
    # plt.ylabel('Susceptibility')
    plt.savefig('cmp.pdf')


if __name__ == '__main__':

    mesh = FDMesh(nx=200, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)

    relax_system(mesh)

    apply_field1(mesh)
    deal_plot()
