import os
from magpy.meshes import CuboidMesh, HexagonalMesh
from magpy.field import scalar_field, vector_field
from magpy.vtk import VTK

MODULE_DIR = os.path.realpath(os.path.dirname(__file__))
OUTPUT_DIR = os.path.join(MODULE_DIR, "temp")


def test_save_scalar_field():
    mesh = CuboidMesh(4, 3, 2, 4, 3, 2)
    s = scalar_field(mesh, lambda r: r[0] + r[1] + r[2])
    vtk = VTK(mesh, directory=OUTPUT_DIR, filename="save_scalar")
    vtk.save_scalar(s, name="s")


def test_save_vector_field():
    mesh = CuboidMesh(4, 3, 2, 4, 3, 2)
    s = vector_field(mesh, lambda r: (r[0], r[1], r[2]))
    vtk = VTK(mesh, directory=OUTPUT_DIR, filename="save_vector")
    vtk.save_vector(s, name="s")


def test_save_scalar_field_hexagonal_mesh():
    mesh = HexagonalMesh(1, 3, 3)
    s = scalar_field(mesh, lambda r: r[0] + r[1])
    vtk = VTK(mesh, directory=OUTPUT_DIR, filename="scalar_hexagonal")
    vtk.save_scalar(s, name="s")
