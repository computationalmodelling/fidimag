import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import FDMesh, UniformExchange, Sim


def init_m(pos):
    x, y, z = pos

    k =  0.1

    nx = k*(x-0.5)

    return (0, np.sin(nx), np.cos(nx))

def test_init():
    mesh = FDMesh(nx=100, ny=1, nz=1)
    sim = Sim(mesh)
    sim.set_m(init_m)

    expected=np.array([0,0,1, 0, np.sin(0.1), np.cos(0.1)])

    assert max(abs(sim.spin[:6]-expected))<1e-15

def test_exch_1d(do_plot=False):
    mesh = FDMesh(nx=100, ny=1, nz=1)
    sim = Sim(mesh)
    sim.set_m(init_m)

    mu0 = 4*np.pi*1e-7
    sim.Ms = 1.0/mu0
    
    exch = UniformExchange(1)
    sim.add(exch)

    field = exch.compute_field()
    field.shape = (-1, 3)

    assert max(abs(field[:,0])) == 0

    xs = np.linspace(0, 99, 100)
    epy = -0.02*np.sin(0.1*xs)
    epz = -0.02*np.cos(0.1*xs)

    assert max(abs(epy[1:-1] - field[1:-1,1]))<3e-5

    if do_plot:
        plt.plot(xs, field[:,1], "-.", label="my", color='DarkGreen')
        plt.plot(xs, field[:,2], "-.", label="mz", color='DarkGreen')
        plt.plot(xs, epy, "--", label="analytical", color='b')
        plt.plot(xs, epz, "--", color='r')
        plt.xlabel("xs")
        plt.ylabel("field")
        plt.legend()
        plt.savefig("exchange_field.pdf") 

if __name__ == '__main__':
    test_init()
    test_exch_1d(do_plot=True)
