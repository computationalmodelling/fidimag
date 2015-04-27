import numpy as np
import scipy.integrate as integrate
from pccp.pc import Constant

const = Constant()

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


"""
compute the figure 2 in the paper at dx.doi.org/10.1063/1.2169472
Note: it seems that the volume given in the paper is not correct, 
since v=28e-27 generates the correct figures. 
"""


def compute_energy(mx):
    Ms = 1.42e6
    a, b, c = 4, 2, 2
    v = 4.0 / 3 * np.pi * a * b * c * 1e-27

    v = 28e-27

    alpha = 0.005
    T = 300
    nu2 = alpha * 2 * const.k_B * T / (const.mu_0 * Ms**2 * v)
    mu = 2 * alpha / nu2

    Dx = 0.0946
    Dy = 0.4132

    fun = lambda t: np.exp(-mu / 2.0 *
                           (Dx * np.cos(t)**2 + Dy * np.sin(t)**2)) * np.sin(t)
    res = integrate.quad(fun, 0, np.pi)
    Z = res[0]

    print 'mu=', mu, 'Z=', Z

    gl = (Dx * mx**2 + Dy * (1 - mx**2)) / 2.0

    p = 1 / Z * np.exp(-mu * gl)

    return p


def plot_energy(name='out.pdf'):
    fig = plt.figure()
    mxs = np.linspace(-1, 1, 21)
    energy = []
    for x in mxs:
        en = compute_energy(x)
        energy.append(en)

    plt.plot(mxs, energy, '.-')
    plt.xlabel('m_x')
    plt.ylabel('Energy')
    #plt.ylim([0, 2])
    fig.savefig(name)


if __name__ == '__main__':
    print compute_energy(0)
    plot_energy()
