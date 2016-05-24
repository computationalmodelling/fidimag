import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, DMI, Zeeman

mu0 = 4 * np.pi * 1e-7


def init_m(pos):

    x, y = pos[0] - 100, pos[1] - 100

    if x**2 + y**2 > 20**2:
        return (0, 0, 1)
    else:
        return (0, 0, -1)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-6, atol=1e-6)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_precession = False

    sim.set_m(init_m)

    exch = UniformExchange(A=1.3e-11)
    sim.add(exch)

    dmi = DMI(D=-4e-3)
    sim.add(dmi)

    zeeman = Zeeman((0, 0, 4e5))
    sim.add(zeeman, save_field=True)

    sim.relax(dt=1e-13, stopping_dmdt=1e-2,
              save_m_steps=None, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)


def excite_system(mesh):
    sim = Sim(mesh, name='dyn', driver='llg_stt')
    sim.set_tols(rtol=1e-8, atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5

    sim.set_m(np.load('m0.npy'))

    exch = UniformExchange(A=1.3e-11)
    sim.add(exch)
    dmi = DMI(D=-4e-3)
    sim.add(dmi)
    zeeman = Zeeman((0, 0, 4e5))
    sim.add(zeeman, save_field=True)

    sim.jx = -5e12
    sim.beta = 0

    ts = np.linspace(0, 0.5e-9, 101)
    for t in ts:
        print 'time', t
        sim.run_until(t)
        sim.save_vtk()
        #sim.save_m()

if __name__ == '__main__':

    mesh = CuboidMesh(nx=101, ny=101, nz=1, dx=2.0, dy=2.0, dz=2.0, unit_length=1e-9, periodicity=(True, True, False))

    relax_system(mesh)

    excite_system(mesh)
