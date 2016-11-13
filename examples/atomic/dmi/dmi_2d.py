import numpy as np
from fidimag.atomistic import Sim, DMI, UniformExchange
from fidimag.common import CuboidMesh

import time


def init_m(pos):

    x, y, z = pos
    if x < 50:
        return (0, 0, 1)
    elif x > 50:
        return (0, 0, -1)
    else:
        return (0, 1, 0)


def relax_system(mesh):

    sim = Sim(mesh, name='dmi_2d')
    sim.driver.alpha = 0.1
    sim.driver.gamma=1.76e11
    sim.mu_s = 1e-22

    J = 1e-20
    exch = UniformExchange(J)
    sim.add(exch)

    dmi = DMI(0.1 * J)
    sim.add(dmi)

    sim.set_m(init_m)    

    ts = np.linspace(0, 5e-10, 101)
    for t in ts:
        print(t)
        sim.run_until(t)
        #sim.save_vtk()

    return sim.spin

if __name__ == '__main__':

    mesh = CuboidMesh(
        nx=1000, ny=1000, nz=1, dx=0.5, dy=0.5, dz=0.5, unit_length=1e-10)

    m0 = relax_system(mesh)
    print 'relax system done'
