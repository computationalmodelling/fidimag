from __future__ import print_function
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from fidimag.common import CuboidMesh
from fidimag.micro import Sim
from fidimag.micro import Zeeman
from fidimag.micro import UniaxialAnisotropy, ExchangeRKKY
import numpy as np


def init_m(pos):
    x, y, z = pos
    if z<1:
        return (0,0.01,-1)
    else:
        return (0.01, 0, 1)


def spatial_Ms(pos):
    x, y, z = pos
    if z<2 and z>1:
        return 0
    else:
        return 8.6e5


def test_excahnge_rkky(do_plot=False):

    mesh = CuboidMesh(nx=1, ny=1, nz=3, dz=1.0, unit_length=1e-9)

    sim = Sim(mesh, name='spin')
    sim.driver.alpha = 0.001
    sim.driver.gamma = 2.211e5
    sim.Ms = spatial_Ms

    sim.set_m(init_m)

    anis = UniaxialAnisotropy(5e4, axis=(0,0,1))
    sim.add(anis)

    exch = ExchangeRKKY(sigma=-1e-4)
    sim.add(exch)

    ts = np.linspace(0, 3e-9, 1001)

    mx = []
    
    mxt = []

    real_ts = []
    for t in ts:
        sim.driver.run_until(t)
        real_ts.append(sim.driver.t)
        print(sim.driver.t, abs(sim.spin_length()[0] - 1))
        mx.append(sim.spin[0])
        mxt.append(sim.spin[6])


    id = np.argmax(abs(np.fft.fft(mx)))
    freqs = np.fft.fftfreq(len(mx), ts[1]-ts[0])

    mu0 = 4*np.pi*1e-7
    Ms = 8.6e5
    K =  2*5e4/(mu0*Ms)
    sigma = 1e-4/1e-9/(mu0*Ms)
    expected = 2.211e5*(K*(K+2*sigma))**0.5/(2*np.pi)

    if do_plot:
        ts_ns = np.array(real_ts) * 1e9
        plt.plot(ts_ns, mx, label="mx", color='darkslateblue')
        plt.plot(ts_ns, mxt, label="mxt", color='m')

        plt.xlabel("time (ns)")
        plt.ylabel("m")
        plt.legend()
        plt.savefig("test_rkky.pdf")

    print((freqs[id]-expected)/expected)
    print((freqs[id]-expected)/expected<0.01)



if __name__ == '__main__':
    test_excahnge_rkky(do_plot=True)
