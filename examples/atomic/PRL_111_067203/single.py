import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fidimag.atomistic import Sim, FDMesh, DMI, UniformExchange, Zeeman
from fidimag.common import Constant


const = Constant()

"""
If we only consider the exchange, dmi and external field, 
we don't have to consider the lattice constant a.

The typical period is (in the paper it should be h rather than h_bar)
    
    T=2*pi/(gamma*J/mu_s) = 2*pi*h_bar*S/J  = h*S/J

so T = 9.6e-13 s and 0.05*T = 4.8e-14s for S=1 and J=50*k_B

"""


def init_m(pos):
    x, y, z = pos

    x0, y0, r = 70, 25, 10

    m1 = (0.05, 0.01, -1)
    m2 = (0, 0, 1)

    if (x - x0)**2 + (y - y0)**2 < r**2:
        return m1
    else:
        return m2


def relax_system(mesh):

    sim = Sim(mesh, name='relax')
    sim.set_default_options(gamma=const.gamma)
    sim.alpha = 0.5
    sim.mu_s = const.mu_s_1

    sim.set_m(init_m)

    J = 50.0 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.5 * J
    dmi = DMI(D)
    sim.add(dmi)

    Hz = 0.2 * J / const.mu_s_1
    zeeman = Zeeman([0, 0, Hz])
    sim.add(zeeman)

    ONE_DEGREE_PER_NS = 17453292.52

    sim.relax(dt=1e-13, stopping_dmdt=0.01 * ONE_DEGREE_PER_NS,
              max_steps=1000, save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def temperature_gradient(pos):

    x = pos[0]

    return x / 150.0 * 0.25 * 50


def excite_system(mesh):

    sim = Sim(mesh, name='dyn', driver='sllg')
    sim.set_options(dt=2e-14, gamma=const.gamma, k_B=const.k_B)
    sim.alpha = 0.1
    sim.mu_s = const.mu_s_1
    sim.T = temperature_gradient

    sim.set_m(np.load("m0.npy"))

    J = 50.0 * const.k_B
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.5 * J
    dmi = DMI(D)
    sim.add(dmi)

    Hz = 0.2 * J / const.mu_s_1
    zeeman = Zeeman([0, 0, Hz])
    sim.add(zeeman)

    dt = 2e-14 * 50  # 1e-12
    ts = np.linspace(0, 1000 * dt, 501)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        sim.save_m()
        print 'sim t=%g' % t


if __name__ == '__main__':
    mesh = FDMesh(nx=150, ny=50, nz=1,  pbc='2d')
    relax_system(mesh)
    excite_system(mesh)
