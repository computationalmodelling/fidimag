import numpy as np
from pccp.pc import *
import time


def init_m(pos):
    x, y, z = pos
    if x < 10:
        return (1, 0, 0)
    elif x > 10:
        return (-1, 0, 0)
    else:
        return (0, 1, 0)


def pin_fun(t, mesh, spin):
    n = mesh.nxyz
    h0 = 0.1
    omega = 3

    t1 = h0 * np.cos(omega * t)
    t2 = h0 * np.sin(omega * t)
    t0 = np.sqrt(1 - t1 * t1 - t2 * t2)

    if t >= 1000 and t < 500:
        t0 = 1
        t1 = 0
        t2 = 0

    spin[0] = t0
    spin[n] = t1
    spin[n + n] = t2


def relax_system(mesh):
    sim = Sim(mesh, name='relax')
    sim.alpha = 0.5
    sim.gamma = 5
    sim.c = 2

    exch = UniformExchange(1)
    sim.add(exch)

    anis = Anisotropy(0.3)
    sim.add(anis)

    sim.set_m(init_m)

    ts = np.linspace(0, 50, 101)
    for t in ts:
        sim.run_until(t)

    # sim.save_vtk()
    return sim.spin


def spin_wave(mesh, m0, H0=10):
    sim = Sim(mesh, name='sw')
    sim.alpha = 0.01
    sim.gamma = 5
    sim.c = 2
    sim.pin_fun = pin_fun

    exch = UniformExchange(1)
    sim.add(exch)

    anis = Anisotropy(0.3)
    sim.add(anis)

    # zeeman=Zeeman(H0,(-1,0,0))
    # sim.add(zeeman)

    sim.set_m(m0)
    # vs=VisualSpin(sim)
    # vs.init()

    ts = np.linspace(0, 500, 501)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        # vs.update()
        # time.sleep(0.001)


if __name__ == '__main__':

    mesh = FDMesh(nx=20)

    m0 = relax_system(mesh)
    print 'relax system done'
    # spin_wave(mesh,m0)
