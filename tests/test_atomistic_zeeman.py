from fidimag.atomistic import Zeeman
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
import numpy as np


def varying_field(pos):
    return (1.2 * pos[0], 2.3 * pos[1], 0)


def test_zeeman():
    """
    Test the x and y component of the zeeman field
    for the 2nd spin

        6  7  8  9  10 11
        0  1  2  3  4  5
              ^
    """

    mesh = CuboidMesh(nx=5, ny=2, nz=1)

    sim = Sim(mesh)
    sim.set_m((1, 0, 0))

    zeeman = Zeeman(varying_field)
    sim.add(zeeman)

    field = zeeman.compute_field()

    assert field[2 * 3] == 1.2 * (2 + 0.5)
    assert field[2 * 3 + 1] == 2.3 * 0.5


if __name__ == '__main__':
    test_zeeman()
