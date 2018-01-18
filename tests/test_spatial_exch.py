import numpy as np
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import Exchange


def init_m(pos):
    x, y, z = pos
    return (x - 0.5, y - 0.5, z - 0.5)


def spatial_J(pos):
    x, y, z = pos
    if x < 5:
        return (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)
    elif x < 6:
        return (-1.0, 0.3, -1.0, -1.0, -1.0, -1.0)
    elif x > 7:
        return (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    elif x > 6:
        return (0.3, 1.0, 1.0, 1.0, 1.0, 1.0)


def test_exch_1d_spatial():
    """
    Test the x component of the exchange field
    in a 1D mesh, with the spin ordering:

    0 1 2 3 4 5

    """
    mesh = CuboidMesh(nx=12, ny=1, nz=1)
    sim = Sim(mesh)
    exch = Exchange(spatial_J)
    sim.add(exch)

    sim.set_m(init_m, normalise=False)

    field = exch.compute_field()

    assert exch._J[3, 3] == -1.0
    assert exch._J[5, 1] == 0.3
    assert exch._J[6, 0] == 0.3
    assert exch._J[8, 5] == 1.0


if __name__ == '__main__':
    test_exch_1d_spatial()
