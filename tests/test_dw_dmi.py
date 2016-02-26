import numpy as np

from fidimag.common import CuboidMesh
from fidimag.micro import Sim
from fidimag.micro import UniformExchange
from fidimag.micro import UniaxialAnisotropy
from fidimag.micro import DMI

import matplotlib.pyplot as plt

mesh = CuboidMesh(dx=2, nx=150, x0=-150, unit_length=1e-9)


def m_init_dw(pos):

    x = pos[0]

    if x < -10:
        return (1, 0, 0)
    elif x > 10:
        return (-1, 0, 0)
    else:
        return (0, 1, 0)


def analytical(xs, A=1.3e-11, D=4e-4, K=8e4):

    delta = np.sqrt(A / (K - D * D / (4 * A))) * 1e9

    phi = D / (2 * A) * xs * 1e-9

    mx = - np.tanh(xs / delta)
    my = 1.0 / np.cosh(xs / delta) * np.cos(phi)
    mz = 1.0 / np.cosh(xs / delta) * np.sin(phi)
    return mx, my, mz


def save_plot(mxyz, mx, my, mz):
    fig = plt.figure()
    mxyz.shape = (-1, 3)
    xs = np.array([p[0] for p in mesh.pos])
    plt.plot(xs, mxyz[:, 0], '.', label='mx')
    plt.plot(xs, mxyz[:, 1], '.', label='my')
    plt.plot(xs, mxyz[:, 2], '.', label='mz')

    plt.plot(xs, mx, '-')
    plt.plot(xs, my, '-')
    plt.plot(xs, mz, '-')

    plt.ylabel('mxyz')
    plt.xlabel('x (nm)')
    plt.legend(loc=1)
    plt.grid()
    plt.xlim([-150, 150])
    plt.ylim([-1.2, 1.2])
    fig.savefig('dw_dmi.pdf')


def test_dw_dmi(mesh=mesh, do_plot=False):

    Ms = 8.0e5
    sim = Sim(mesh, name='relax')

    sim.set_m(m_init_dw)

    sim.set_tols(rtol=1e-8, atol=1e-12)
    sim.Ms = Ms
    sim.alpha = 0.5
    sim.do_precession = False

    A = 1.3e-11
    D = 4e-4
    Kx = 8e4
    Kp = -6e5

    sim.add(UniformExchange(A))
    sim.add(DMI(D))
    sim.add(UniaxialAnisotropy(Kx, axis=[1, 0, 0], name='Kx'))

    sim.relax(stopping_dmdt=0.01)

    xs = np.array([p[0] for p in mesh.coordinates])
    mx, my, mz = analytical(xs, A=A, D=D, K=Kx)
    mxyz = sim.spin.copy()
    mxyz = mxyz.reshape(-1, 3)

    assert max(abs(mxyz[:, 0] - mx)) < 0.002
    assert max(abs(mxyz[:, 1] - my)) < 0.002
    assert max(abs(mxyz[:, 2] - mz)) < 0.0006

    if do_plot:

        save_plot(mxyz, mx, my, mz)

if __name__ == '__main__':

    test_dw_dmi(do_plot=True)
