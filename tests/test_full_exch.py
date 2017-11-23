from __future__ import print_function
import numpy as np
from fidimag.atomistic import Sim
from fidimag.atomistic import HexagonalMesh
from fidimag.atomistic import Exchange


def init_m(pos, a):
    """
    This function makes the sx, sy, sz components of a hex lattice
    to be the (i, j) index of their position
    """
    x, y, z = pos
    base = a * 0.5 / np.cos(np.pi / 6.)
    dy = a * np.sin(np.pi / 3.)

    new_y = (y - base) / dy
    new_x = x - a * 0.5
    if new_y % 2 == 1:
        new_x += 0.5

    return (int(new_x), int(new_y), z)


def remove_negatives(a):
    """
    Return *a* discarding elements with -1
    """
    return a[a != -1]


def test_full_exch_hex_3_shells_lattice_pos():
    """
    Test the x component of the exchange field when using 3 shells of
    neighbours in a 9 x 9 hexagonal lattice with square arrangement

    We set the s_x and s_y spin components as the (i, j) index of the
    lattice positions:
    """
    a = 1
    shells = 3
    mesh = HexagonalMesh(a * 0.5, nx=9, ny=9,
                         shells=shells, alignment='square')
    sim = Sim(mesh)
    # Set the indexes (i, j) as their s_x, s_y position
    sim.set_m(lambda r: init_m(r, a), normalise=False)

    Js = np.ones(shells)
    exch = Exchange(Js)
    sim.add(exch)

    field = exch.compute_field()
    assert field[3 * 0] == (1 + 1 + 0) + (2 + 0) + (2 + 1)  # f_x 1st spin
    assert field[3 * 11] == ((3 + 1 + 2 + 1 + 1 + 2) +
                             (3 + 0 + 2 + 0 + 0 + 3) +
                             (4 + 0 + 3 + 0 + 1 + 0))


def test_full_exch_hex_9_shells():
    """
    Test the x component of the exchange field when using 9 shells of
    neighbours in a 9 X 9 hexagonal lattice with square and diagonal
    arrangement

    We set J=1 for every shell and set the s_x component of the spins as the
    lattice site number:
    [0, 1, 2, 3, ... 80]

    Since we set the s_x components as the lattice position indexes, the x
    component of the field is just the sum of the indexes of the neighbours
    (assuming the neighbours indexing is correct) because we set the exchange
    constants as 1
    """

    for arrang in ['square', 'diagonal']:
        a = 1
        shells = 9
        mesh = HexagonalMesh(a * 0.5, nx=9, ny=9,
                             shells=shells, alignment=arrang)
        sim = Sim(mesh)
        # Set s_x as the lattice site number
        sim.spin.reshape(-1, 3)[:, 0] = np.arange(len(sim.spin.reshape(-1, 3)[:, 0]))

        Js = np.ones(shells)
        exch = Exchange(Js)
        sim.add(exch)

        field = exch.compute_field()

        # Test the x component for every lattice site, summing the neighbours
        # index and removing the -1 elements (their contribution is zero)
        for i in range(sim.mesh.n):
            assert field[3 * i] == np.sum(remove_negatives(sim.mesh.neighbours[i]))


def test_full_exch_hex_9_shells_J_rings():
    """
    Test the x component of the exchange field when using 9 shells of
    neighbours in a 11 X 11 hexagonal lattice with square and diagonal
    arrangement

    We set J=1,2,3,.. for every shell and set the s_x component of the spins as
    the lattice site number:
    [0, 1, 2, 3, ... 120]

    Since we set the s_x components as the lattice position indexes, the x
    component of the field is the sum of the indexes of the neighbours
    (assuming the neighbours indexing is correct) multiplied by J[i] where i is
    the shell (1, 2, ...9), i.e. the 1st shell of ngbs is multiplied by 1,
    the 2nd shell by 2, the 3rd shell by 3, and so on
    """

    for arrang in ['square', 'diagonal']:
        a = 1
        shells = 9
        mesh = HexagonalMesh(a * 0.5, nx=11, ny=11,
                             shells=shells, alignment=arrang)
        sim = Sim(mesh)
        # Set s_x as the lattice site number
        sim.spin.reshape(-1, 3)[:, 0] = np.arange(len(sim.spin.reshape(-1, 3)[:, 0]))

        # Exchange constants according to the shell
        Js = np.arange(1, 10)
        exch = Exchange(Js)
        sim.add(exch)

        field = exch.compute_field()

        # We only test for the 60th lattice site
        ngbs_60 = mesh.neighbours[60]
        sum_ngbs = mesh._sum_ngbs_shell

        # For every shell, find the ngb indexes in that shell and multiply the
        # sum by the corresponding J=1, 2, 3, ...
        sum_rings = 0
        for i in range(1, shells + 1):
            ngbs_range = slice(sum_ngbs[i - 1], sum_ngbs[i])
            print('J = ', Js[i - 1], '  ngbs indexes: ', ngbs_60[ngbs_range])
            sum_rings += Js[i - 1] * np.sum(ngbs_60[ngbs_range])

        assert field[3 * 60] == sum_rings


def test_full_exch_hex_2_shells():
    """
    Test the x component of the exchange field when using 2 shells of
    neighbours, comparing the field manually.
    This is similar than the *test_full_exch_hex_9_shells* function but
    here we do not assume that the neighbours indexes are correct
    We set J=1 for NN and J=2 for NNN
    """
    a = 1
    shells = 2
    mesh = HexagonalMesh(a * 0.5, nx=9, ny=9,
                         shells=shells, alignment='square')
    sim = Sim(mesh)
    sim.spin.reshape(-1, 3)[:, 0] = np.arange(len(sim.spin.reshape(-1, 3)[:, 0]))
    Js = np.array([1., 2.])
    exch = Exchange(Js)
    sim.add(exch)
    field = exch.compute_field()
    assert field[3 * 0] == 1. * (1 + 10 + 9) + 2. * (11 + 18)
    assert field[3 * 11] == (12 + 10 + 20 + 1 + 19 + 2) + 2. * (21 + 29 + 18 + 3)


if __name__ == '__main__':
    test_full_exch_hex_2_shells()
    test_full_exch_hex_9_shells()
    test_full_exch_hex_9_shells_J_rings()
