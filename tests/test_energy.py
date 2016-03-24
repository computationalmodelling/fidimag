from __future__ import print_function
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Anisotropy
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
from fidimag.atomistic import UnitMaterial
from fidimag.common.fileio import DataReader
import numpy as np


class UnitMaterial(object):

    def __init__(self):
        self.a = 1
        self.b = 1
        self.c = 1
        self.J = 1.0
        self.Dx = 0.5
        self.mu_s = 1.0
        self.gamma = 2 * np.pi
        self.alpha = 0.01
        self.unit_length = 1.0


def pin_fun(pos):
    if pos[0] == 0:
        return 1
    else:
        return 0


def relax_system(mesh, Dx=0.005, Dp=0.01):

    mat = UnitMaterial()

    sim = Sim(mesh, name='test_energy')
    print('Created sim')
    sim.set_tols(rtol=1e-10, atol=1e-12)

    sim.alpha = mat.alpha
    sim.gamma = mat.gamma
    sim.pins = pin_fun

    exch = UniformExchange(mat.J)
    sim.add(exch)
    print('Added UniformExchange')

    anis = Anisotropy(Dx, axis=[1, 0, 0], name='Dx')
    sim.add(anis)
    print('Added Anisotropy')

    anis2 = Anisotropy([0, 0, -Dp], name='Dp')
    sim.add(anis2)
    print('Added Anisotropy 2')

    sim.set_m((1, 1, 1))

    T = 100
    ts = np.linspace(0, T, 201)
    for t in ts:
        # sim.save_vtk()
        sim.run_until(t)
        print(('Running -', t))

    # sim.save_vtk()
    np.save('m0.npy', sim.spin)


def save_plot():
    fig = plt.figure()

    data = DataReader('test_energy.txt')
    ts = data['time']

    #plt.plot(ts, data['E_Dp'], label='E_Dp')
    #plt.plot(ts, data['E_Dp'], label='E_Dp')
    #plt.plot(ts, data['E_Dx'], label='E_Dx')
    plt.plot(ts, data['E_exch'], label='E_exch')
    plt.plot(ts, data['E_total'], label='E_total')

    plt.ylabel('Energy (J)')
    plt.xlabel('Time (s)')
    plt.legend(loc=1)
    plt.grid()
    # plt.ylim([-1.2,1.2])
    fig.savefig('energy.pdf')


def test_energy(do_plot=False):
    mesh = CuboidMesh(nx=30)
    print('Made mesh')
    relax_system(mesh)
    print('Done relax_system')

    if do_plot:
        save_plot()

    data = DataReader('test_energy.txt')
    print('Made DataReader')
    energy = data['E_total']

    for i in range(len(energy) - 1):
        assert energy[i] > energy[i + 1]


if __name__ == '__main__':

    test_energy(do_plot=False)
