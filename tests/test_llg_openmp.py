from __future__ import print_function
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from fidimag.common import CuboidMesh
from fidimag.micro import Sim
from fidimag.micro import Zeeman
from fidimag.micro import UniaxialAnisotropy
import numpy as np


def init_m(pos):
    x, y, z = pos
    if (x, y, z) == (1 + 0.5, 2 + 0.5, 3 + 0.5):
        return (1, 2, 3)
    elif z < 1 + 0.5:
        return (0, 0, -1)
    else:
        return (0, 0, 1)


def pin_fun(pos):
    if pos[0] < 1:
        return 1
    else:
        return 0


def test_sim_pin():
    mesh = CuboidMesh(nx=3, ny=2, nz=1)
    sim = Sim(mesh, integrator='sundials_openmp')
    sim.set_m((0, 0.8, 0.6))
    sim.alpha = 0.1
    sim.driver.gamma = 1.0
    sim.pins = pin_fun

    anis = UniaxialAnisotropy(Ku=1, axis=[0, 0, 1], name='Dx')
    sim.add(anis)

    sim.driver.run_until(1.0)
    print(sim.spin)
    assert sim.spin[0] == 0
    assert sim.spin[2] != 0


def test_sim_init_m():
    mesh = CuboidMesh(nx=3, ny=4, nz=5)
    sim = Sim(mesh)
    sim.set_m((0, 1, 0))
    sim.spin.shape = (-1, 3)
    spin_y = sim.spin[:, 1]
    assert(spin_y.any() == 1)


def test_sim_init_m_fun():
    mesh = CuboidMesh(nx=3, ny=4, nz=5)
    sim = Sim(mesh)
    sim.set_m(init_m, normalise=False)

    print(sim.spin.reshape(-1, 3).shape)
    print(sim.mesh.index(1, 2, 3))

    assert(sim.spin_at(1, 2, 3)[0] == 1)
    assert(sim.spin_at(1, 2, 3)[1] == 2)
    assert(sim.spin_at(1, 2, 3)[2] == 3)


def test_m_average():
    mesh = CuboidMesh(nx=3, ny=4, nz=5)
    sim = Sim(mesh)
    sim.set_m((0, 0, 1))
    a = sim.compute_average()
    assert a[2] == 1.0


def single_spin(alpha, gamma, H0, ts):
    """
    compute single spin under the external field H
    """

    precession = gamma / (1 + alpha**2)
    beta = precession * H0 * ts

    mx = np.cos(beta) / np.cosh(alpha * beta)
    my = np.sin(beta) / np.cosh(alpha * beta)
    mz = np.tanh(alpha * beta)

    return mx, my, mz


def test_sim_single_spin(do_plot=False):

    mesh = CuboidMesh(nx=80, ny=3, nz=3)

    sim = Sim(mesh, name='spin', integrator='sundials_openmp')

    alpha = 0.1
    gamma = 2.21e5
    sim.alpha = alpha
    sim.driver.gamma = gamma
    sim.mu_s = 1.0

    sim.set_m((1, 0, 0))

    H0 = 1e5
    sim.add(Zeeman((0, 0, H0)))

    ts = np.linspace(0, 1e-9, 101)

    mx = []
    my = []
    mz = []
    real_ts = []
    for t in ts:
        sim.driver.run_until(t)
        real_ts.append(sim.driver.t)
        print(sim.driver.t, abs(sim.spin_length()[0] - 1))
        mx.append(sim.spin[0])
        my.append(sim.spin[1])
        mz.append(sim.spin[2])

    mz = np.array(mz)
    # print mz
    a_mx, a_my, a_mz = single_spin(alpha, gamma, H0, ts)

    print(sim.driver.stat())

    if do_plot:
        ts_ns = np.array(real_ts) * 1e9
        plt.plot(ts_ns, mx, ".", label="mx", color='DarkGreen')
        plt.plot(ts_ns, my, ".", label="my", color='darkslateblue')
        plt.plot(ts_ns, mz, ".", label="mz", color='m')
        plt.plot(ts_ns, a_mx, "--", label="analytical", color='b')
        plt.plot(ts_ns, a_my, "--",  color='b')
        plt.plot(ts_ns, a_mz, "--",  color='b')
        plt.xlabel("time (ns)")
        plt.ylabel("m")
        plt.title("integrating a macrospin")
        plt.legend()
        plt.savefig("single_spin.pdf")

    print(("Max Deviation = {0}".format(
        np.max(np.abs(mz - a_mz)))))

    assert np.max(np.abs(mz - a_mz)) < 5e-7


if __name__ == '__main__':
    test_sim_single_spin(do_plot=True)
