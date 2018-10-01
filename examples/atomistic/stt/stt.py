import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange, Zeeman, Anisotropy
from fidimag.common import CuboidMesh


def init_m(pos):
    x, y, z = pos
    if x < 250:
        return (1, 0, 0)
    elif x > 270:
        return (-1, 0, 0)
    else:
        return (0, 1, 0)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')
    # sim.set_options(rtol=1e-10,atol=1e-14)
    sim.driver.alpha = 1.0
    sim.driver.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m(init_m)
    # sim.set_m(random_m)
    # sim.set_m(np.load('m_10000.npy'))

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    Kx = Anisotropy(Ku=0.005, axis=(1, 0, 0), name='Kx')
    sim.add(Kx)

    sim.relax(dt=2.0, stopping_dmdt=1e-6, max_steps=1000,
              save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def dynamic(mesh):

    sim = Sim(mesh, name='dyn', driver='slonczewski')
    # sim.set_options(rtol=1e-10,atol=1e-14)
    sim.driver.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m(np.load('m0.npy'))

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    Kx = Anisotropy(Ku=0.005, axis=(1, 0, 0), name='Kx')
    sim.add(Kx)

    sim.p = (0,0,1)

    sim.u0 = 0.03
    sim.driver.alpha = 0.1

    ts = np.linspace(0, 1e3, 101)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        print t


if __name__ == '__main__':
    mesh = CuboidMesh(nx=500, ny=1, nz=1)
    relax_system(mesh)
    dynamic(mesh)
