import numpy as np
from micro import Sim
from micro import FDMesh
from micro import Zeeman, TimeZeeman
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
    mesh = FDMesh(nx=1, ny=1, nz=1)
    sim = Sim(mesh, name='relax')
    sim.set_tols(rtol=1e-10, atol=1e-10)
    sim.alpha = 0.1

    sim.set_m((0, 0, 1))

    sim.add(Zeeman((0, 0, 1e5)))

    w0 = 2 * np.pi * 1e9

    def sin_fun(t):
        return np.sin(w0 * t)

    h0 = 1e3
    theta = np.pi / 20
    hx = h0 * np.sin(theta)
    hz = h0 * np.cos(theta)
    hx = TimeZeeman([hx, 0, hz], sin_fun, name='h')
    sim.add(hx, save_field=True)

    ts = np.linspace(0, 5e-9, 5001)

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
    #plt.plot(ts, mx, '--', label='m_x',dashes=(2.0,2.0))
    #plt.plot(ts, my, '--', label='m_y',dashes=(2.0,2.0))
    #plt.plot(ts, mz, '--', label='m_z',dashes=(2.0,2.0))

    ts = ts[::10]
    mx = data['m_x'][::10]
    my = data['m_y'][::10]
    mz = data['m_z'][::10]
    plt.plot(ts, mx, color='b', label='m_x')
    plt.plot(ts, my, color='g', label='m_y')
    # plt.plot(ts[::2],mz[::2],color='r')
    # plt.xlim([0,1.01])

    l1 = plt.legend(bbox_to_anchor=[0.8, 0.8], shadow=True, frameon=True)
    custom_legend(l1)

    plt.xlabel('Time (ns)')
    plt.ylabel('m')
    plt.tight_layout()
    fig.savefig('m_ts.pdf')


if __name__ == '__main__':

    relax_system()
    plot_all()
