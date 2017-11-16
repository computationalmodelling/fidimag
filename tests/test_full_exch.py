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


def test_full_exch_hex():
    """

    """
    a = 1
    mesh = HexagonalMesh(a * 0.5, nx=9, ny=9, shells=2)
    sim = Sim(mesh)
    sim.set_m(lambda r: init_m(r, a), normalise=False)
    sim.spin.reshape(-1, 3)[:, 0] = np.arange(len(sim.spin.reshape(-1, 3)[:, 0]))

    # print(sim.mesh.coordinates[:])
    # print(sim.spin.reshape(-1, 3)[:])

    Js = [1, 1]
    exch = Exchange(Js)
    sim.add(exch)
    exch.compute_field()

    print(field)
    print(sim.get_interaction('Exchange').field.reshape(-1, 3)[:5])

    # assert field[3 * 40] =


if __name__ == '__main__':
    test_full_exch_hex()
