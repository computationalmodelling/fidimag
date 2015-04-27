import numpy as np
from fidimag.pc import Sim
from fidimag.pc import FDMesh
from fidimag.pc import UniformExchange


def init_m(pos):
    x, y, z = pos
    return (x - 0.5, y - 0.5, z - 0.5)


def test_exch_1d():
    mesh = FDMesh(nx=5, ny=1, nz=1)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    assert field[0] == 1
    assert field[1] == 2
    assert field[2] == 4
    assert field[3] == 6
    assert field[4] == 3
    assert np.max(field[5:]) == 0


def test_exch_1d_pbc():
    mesh = FDMesh(nx=5, ny=1, nz=1, pbc='x')
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    assert field[0] == 1 + 4
    assert field[1] == 2
    assert field[2] == 4
    assert field[3] == 6
    assert field[4] == 3 + 0
    assert np.max(field[5:]) == 0


def test_exch_2d():
    mesh = FDMesh(nx=5, ny=2, nz=1)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    assert field[0] == 1
    assert field[2] == 2 + 1
    assert field[4] == 1 + 2 + 3
    assert field[6] == 2 + 3 + 4
    assert field[8] == 3 + 4
    assert np.max(field[21:]) == 0


def test_exch_2d_pbc2d():
    mesh = FDMesh(nx=3, ny=2, nz=1, pbc='xy')
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()

    expected_x = np.array([3, 3, 4, 4, 5, 5])
    expected_y = np.array([2, 2, 2, 2, 2, 2])

    assert np.max(abs(field[:6] - expected_x)) == 0
    assert np.max(abs(field[6:12] - expected_y)) == 0

    assert np.max(field[12:]) == 0


def test_exch_3d():
    mesh = FDMesh(nx=4, ny=3, nz=2)
    sim = Sim(mesh)
    exch = UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()
    # print field
    assert field[0] == 1
    assert field[6] == 0 + 1 + 2 + 1
    assert field[12] == 1 + 2 + 3 + 2
    assert field[18] == 2 + 3 + 3

    assert field[2] == 1
    assert field[8] == 5
    assert field[14] == 10
    assert field[20] == 11


def test_exch_energy_1d():
    mesh = FDMesh(nx=2, ny=1, nz=1)
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
