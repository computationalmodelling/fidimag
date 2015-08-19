import numpy as np
from fidimag.atomistic import Sim, FDMesh, DMI, UniformExchange


def init_m(pos):
    x, y, z = pos
    if x < 50:
        return (0, 0, 1)
    elif x > 50 - 1:
        return (0, 1, -1)
    else:
        return (0, 1, 0)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')
    sim.alpha = 0.1

    sim.set_m(init_m)

    J = 1
    exch = UniformExchange(J)
    sim.add(exch)

    dmi = DMI(0.05 * J)
    sim.add(dmi)

    ts = np.linspace(0, 1, 11)
    for t in ts:
        print t, sim.spin_length() - 1
        sim.run_until(t)

    sim.save_vtk()

    return sim.spin


if __name__ == '__main__':

    mesh = FDMesh(nx=20, ny=1, nz=1, pbc=None)

    m0 = relax_system(mesh)
    print 'relax system done'
    # spin_wave(mesh,m0)
