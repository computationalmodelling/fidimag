import numpy as np
from field import scalar_field, vector_field
from cuboid_mesh import CuboidMesh


def test_initialise_scalar():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    f = scalar_field(mesh, lambda r: r[0] + r[1] + r[2])
    assert np.allclose(f, np.array((1.5)))


def test_initialise_vector():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    v = vector_field(mesh, lambda r: 2 * r)
    assert np.allclose(v, np.array((1, 1, 1)))


