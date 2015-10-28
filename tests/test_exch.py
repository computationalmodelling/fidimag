import numpy as np
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import UniformExchange


def init_m(pos):
    x, y, z = pos
    return (x - 0.5, y - 0.5, z - 0.5)


def test_exch_1d():
    """
    Test the x component of the exchange field
    in a 1D mesh, with the spin ordering:

    0 1 2 3 4 5

    """
    mesh = CuboidMesh(nx=5, ny=1, nz=1)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()

    assert field[0] == 1
    assert field[1 * 3] == 2
    assert field[2 * 3] == 4
    assert field[3 * 3] == 6
    assert field[4 * 3] == 3

    assert np.max(field[2::3]) == 0
    assert np.max(field[1::3]) == 0


def test_exch_1d_pbc():
    mesh = CuboidMesh(nx=5, ny=1, nz=1, periodicity=(True, False, False))
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    assert field[0] == 1 + 4
    assert field[3] == 2
    assert field[6] == 4
    assert field[9] == 6
    assert field[12] == 3 + 0
    assert np.max(field[2::3]) == 0
    assert np.max(field[1::3]) == 0


def test_exch_2d():
    mesh = CuboidMesh(nx=5, ny=2, nz=1)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()

    assert np.max(field[2::3]) == 0

    assert field[0] == 1
    assert field[3] == 2 + 1
    assert field[6] == 1 + 2 + 3
    assert field[9] == 2 + 3 + 4
    assert field[12] == 3 + 4


def test_exch_2d_pbc2d():
    """
    Test the exchange field components in a 2D mesh with PBCs
    The mesh sites:

            3     4     5    -->    (0,1,0)  (1,1,0)  (2,1,0)
     y ^    0     1     2           (0,0,0)  (1,0,0)  (2,0,0)
       |
       x -->
  
    The expected components are in increasing order along x

    """

    mesh = CuboidMesh(nx=3, ny=2, nz=1, periodicity=(True, True, False))
    print mesh.neighbours
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()

    expected_x = np.array([3, 4, 5, 3, 4, 5])
    expected_y = np.array([2, 2, 2, 2, 2, 2])

    # Since the field ordering is now: fx1 fy1 fz1 fx2 ...
    # We extract the x components jumping in steps of 3
    assert np.max(abs(field[::3] - expected_x)) == 0
    # For the y component is similar, now we start at the 1th
    # entry and jump in steps of 3
    assert np.max(abs(field[1::3] - expected_y)) == 0
    # Similar fot he z component
    assert np.max(field[2::3]) == 0


def test_exch_3d():
    """
    Test the exchange field of the spins in this 3D mesh:

    bottom layer:
    8  9  10  11
    4  5  6   7       x 2
    0  1  2   3

    The assertions are the mx component
    of the: 0, 1, 2, .. 7 spins

    Remember the new new ordering: fx1, fy1, fz1, fx2, ...

    """
    mesh = CuboidMesh(nx=4, ny=3, nz=2)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    # print field
    assert field[0] == 1
    assert field[3] == 0 + 1 + 2 + 1
    assert field[6] == 1 + 2 + 3 + 2
    assert field[9] == 2 + 3 + 3

    assert field[4 * 3] == 1
    assert field[5 * 3] == 5
    assert field[6 * 3] == 10
    assert field[7 * 3] == 11


def test_exch_energy_1d():
    mesh = CuboidMesh(nx=2, ny=1, nz=1)
    sim = Sim(mesh)
    exch = UniformExchange(1.23)
    sim.add(exch)

    sim.set_m((0, 0, 1))

    energy = exch.compute_energy()
    assert energy == -1.23


if __name__ == '__main__':
    # test_exch_1d()
    # test_exch_1d_pbc()
    # test_exch_2d()
    test_exch_2d_pbc2d()
    # test_exch_3d()
    # test_exch_energy_1d()
