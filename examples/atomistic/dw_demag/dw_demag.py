import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fidimag.atomistic import Sim
from fidimag.atomistic import DMI, Demag, UniformExchange
from fidimag.common import CuboidMesh


def init_m(pos):
    x, y, z = pos
    if x < 140 * 0.5:
        return (1, 0, 0.1)
    elif x > 160 * 0.5:
        return (-1, 0, 0.1)
    else:
        return (0, 1, 0.1)


def relax_system(mesh):

    sim = Sim(mesh, name='relax')
    sim.mu_s=1e-23
    sim.driver.gamma=1.76e11
    sim.driver.alpha = 1.0
    J = 1e-22
    exch = UniformExchange(J)
    sim.add(exch)
    demag = Demag()
    sim.add(demag)
    sim.set_m(init_m)

    ts = np.linspace(0, 5e-10, 101)
    for t in ts:
        sim.driver.run_until(t)
        sim.save_vtk()
    np.save('m0.npy', sim.spin)


def save_plot():
    fig = plt.figure()
    data = np.load('m0.npy')
    data.shape = (3, -1)
    print(data)

    plt.plot(data[0], '-', label='Sx')
    plt.plot(data[1], '-', label='Sy')
    plt.plot(data[2], '-', label='Sz')

    plt.ylim([-1.2, 1.2])
    plt.ylabel('S')
    plt.xlabel('x/a')
    plt.legend(loc=1)
    plt.grid()
    fig.savefig('mxyz.pdf')

if __name__ == '__main__':

    mesh = CuboidMesh(nx=300, dx=0.5, dy=1, dz=1, unit_length=1e-9)

    relax_system(mesh)
    print('relax system done')
    save_plot()
    # spin_wave(mesh,m0)
