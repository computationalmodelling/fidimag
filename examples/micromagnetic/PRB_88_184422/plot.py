import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim
from micro import FDMesh
from micro import Zeeman
from fidimag.common.fileio import DataReader

from util.omf import OMF2

mesh = FDMesh(nx=100, dx=1, unit_length=1e-9)


def plot_all():

    font = {'family': 'serif',
            'weight': 'normal',
            'size': 12,
            }

    plt.rc('font', **font)

    Ms = 8.6e5
    omf = OMF2('oommf/dmi-Oxs_TimeDriver-Magnetization-00-0000963.omf')
    mx = omf.get_all_mag(comp='x') / Ms
    my = omf.get_all_mag(comp='y') / Ms
    mz = omf.get_all_mag(comp='z') / Ms
    xs = np.linspace(0.5, 100 - 0.5, 100)

    fig = plt.figure(figsize=(5, 4))
    plt.plot(xs, mx, '--', label='m_x', dashes=(2.0, 2.0))
    plt.plot(xs, my, '--', label='m_y', dashes=(2.0, 2.0))
    plt.plot(xs, mz, '--', label='m_z', dashes=(2.0, 2.0))

    l1 = plt.legend(bbox_to_anchor=[0.8, 0.8], shadow=True, frameon=True)
    # custom_legend(l1)

    plt.xlabel('xs (nm)')
    plt.ylabel('m')
    plt.title('OOMMF')
    plt.tight_layout()
    fig.savefig('m_oommf.pdf')


if __name__ == '__main__':

    # relax_system()
    plot_all()
