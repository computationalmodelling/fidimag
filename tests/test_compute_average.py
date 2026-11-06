import numpy as np

from fidimag.common import CuboidMesh
from fidimag.micro import Sim


def test_the_average_is_over_the_cells_with_material():
    """`compute_average` divided by the number of cells, so a mesh with empty cells gave an average that was too small"""
    mesh = CuboidMesh(nx=4, ny=1, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh)
    sim.set_Ms(lambda r: 0.0 if r[0] > 4 else 8e5)     # the last two cells are empty
    sim.set_m((0, 0, 1))
    assert np.allclose(sim.compute_average(), (0, 0, 1))


def test_the_average_of_a_full_mesh_is_unchanged():
    mesh = CuboidMesh(nx=3, ny=2, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh)
    sim.set_Ms(8e5)
    sim.set_m(lambda r: (1, 0, 0) if r[0] < 3 else (-1, 0, 0))
    assert np.allclose(sim.compute_average(), (-1 / 3, 0, 0))     # (two cells at +1, four at -1)
