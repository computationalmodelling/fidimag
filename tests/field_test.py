import numpy as np
from magpy.field import scalar_field, vector_field
from magpy.meshes.cuboid_mesh import CuboidMesh
from . import very_close


def test_initialise_scalar():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    f = scalar_field(mesh, lambda r: r[0] + r[1] + r[2])
    assert very_close(f, np.array((1.5)))


def test_initialise_vector():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    v = vector_field(mesh, lambda r: 2 * r)
    assert very_close(v, np.array((1, 1, 1)))


