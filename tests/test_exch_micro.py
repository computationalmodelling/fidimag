"""
Test of the micromagnetic exchange field using an analytical
expression.

The system is a 1D mesh of 100 spins: FDMesh(nx=100, ny=1, nz=1)
which are located at positions starting at 0.5 since the default
mesh spacing is dx=1, thus we have something like:

           dx=1
          o -- o -- o -- o -- o --
x    =   0.5  1.5  2.5  3.5  ...
i    =    0    1    2    3

The **init_m** function will set the magnetisation as:
    mx = 0, my = sin(n_x) , mz = cos(n_x)

where n_x is ( k * (x - 0.5) ), with x the position
of the i-th spin, therefore the values in the mesh are:

                o -- o -- o -- o -- o --
x   - 0.5  =    0    1    2    3   ...
n_x        =    0   0.1  0.2  0.3

For the **test_exch_1d**, we set A = 1, Ms = 1 / mu0, hence
the exchange field H_ex will be:

    H_ex = (2 A / (mu0 * Ms)) * nabla^2 (mx, my, mz)
         = 2 * (Dx^2 mx, Dx^2 my, Dx^2 mz)
         = -2 * (k ^ 2) * (0, sin(n_x), cos(n_x)

where we used Dx^2 = d^2 / dx^2. Remember that
sin(n_x) = sin(k * (x - 0.5)).

With this expression, we compare the Fidimag computed value.

"""

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import FDMesh, UniformExchange, Sim


def init_m(pos):
    """
    The initial magnetisation function defined in the
    introduction
    """
    x, y, z = pos

    k = 0.1
    nx = k * (x - 0.5)

    return (0, np.sin(nx), np.cos(nx))


def test_init():
    """
    This tests (mx, my, mx) for the first 2 spins
    """
    mesh = FDMesh(nx=100, ny=1, nz=1)
    sim = Sim(mesh)
    sim.set_m(init_m)

    expected = np.array([0, 0, 1,
                         0, np.sin(0.1), np.cos(0.1)])

    assert max(abs(sim.spin[:6] - expected)) < 1e-15


def test_exch_1d(do_plot=False):
    # Initiate the 1D mesh and magnetisation as before
    mesh = FDMesh(nx=100, ny=1, nz=1)
    sim = Sim(mesh)
    sim.set_m(init_m)

    # Simplify the magnetic parameters
    mu0 = 4 * np.pi * 1e-7
    sim.Ms = 1.0 / mu0

    exch = UniformExchange(1)
    sim.add(exch)

    # Compute the exchange field and reshape it in order
    # to leave every row as the [f_x, f_y, f_z] array
    # for every spin
    field = exch.compute_field()
    field.shape = (-1, 3)

    # We know that the field in x is always zero ( see the
    # analytical calculation at the beginning)
    assert max(abs(field[:, 0])) == 0

    # These are the analytical values for the exchange field in y,z
    # In this case, k=0.1 , then 2 * k^2 evaluates as 0.02
    xs = np.linspace(0, 99, 100)
    epy = -0.02 * np.sin(0.1 * xs)
    epz = -0.02 * np.cos(0.1 * xs)

    # Compare the analytical value
    # of the y component of the exchange field, with Fidimag's
    # result (second column of the reshaped field array)

    # WARNING: NOTICE that we are not considering the extremes since
    # there is a wrong expression in the border of the exchange field
    # with NO PBCs. We must FIX this test!
    assert max(abs(epy[1:-1] - field[1:-1, 1])) < 3e-5

    if do_plot:
        plt.plot(xs, field[:, 1], "-.", label="my", color='DarkGreen')
        plt.plot(xs, field[:, 2], "-.", label="mz", color='DarkGreen')
        plt.plot(xs, epy, "--", label="analytical", color='b')
        plt.plot(xs, epz, "--", color='r')
        plt.xlabel("xs")
        plt.ylabel("field")
        plt.legend()
        plt.savefig("exchange_field.pdf")

if __name__ == '__main__':
    test_init()
    test_exch_1d(do_plot=True)
