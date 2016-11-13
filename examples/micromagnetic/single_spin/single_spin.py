import numpy as np
from micro import Sim
from common import CuboidMesh
from micro import Zeeman
from fidimag.common.fileio import DataReader
import matplotlib.pyplot as plt


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


def relax_system():
    mesh = CuboidMesh(nx=1, ny=1, nz=1)
    sim = Sim(mesh, name='relax')
    sim.driver.set_tols(rtol=1e-10, atol=1e-10)
    sim.driver.alpha = 0.5

    sim.set_m((1.0, 0, 0))

    sim.add(Zeeman((0, 0, 1e5)))

    ts = np.linspace(0, 1e-9, 1001)

    for t in ts:
        sim.run_until(t)


def custom_legend(legend):
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize(11)

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width


def plot_all():

    font = {'family': 'serif',
            'weight': 'normal',
            'size': 12,
            }

    plt.rc('font', **font)

    data = DataReader('relax.txt')
    ts = data['time']

    mx, my, mz = single_spin(0.5, 2.21e5, 1e5, ts)

    ts = ts * 1e9
    fig = plt.figure(figsize=(5, 4))
    plt.plot(ts, mx, '--', label='m_x', dashes=(2.0, 2.0))
    plt.plot(ts, my, '--', label='m_y', dashes=(2.0, 2.0))
    plt.plot(ts, mz, '--', label='m_z', dashes=(2.0, 2.0))

    ts = ts[::10]
    mx = data['m_x'][::10]
    my = data['m_y'][::10]
    mz = data['m_z'][::10]
    plt.plot(ts, mx, '.', color='b')
    plt.plot(ts, my, '.', color='g')
    plt.plot(ts[::2], mz[::2], '.', color='r')
    plt.xlim([0, 1.01])

    l1 = plt.legend(bbox_to_anchor=[0.8, 0.8], shadow=True, frameon=True)
    custom_legend(l1)

    plt.xlabel('Time (ns)')
    plt.ylabel('m')
    plt.tight_layout()
    fig.savefig('m_ts.pdf')


if __name__ == '__main__':

    relax_system()
    plot_all()
