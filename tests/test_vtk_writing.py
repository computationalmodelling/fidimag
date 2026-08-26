import os
import base64
import struct
import xml.etree.ElementTree as ET

import numpy as np

from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.common.field import scalar_field, vector_field
from fidimag.common.vtk import VTK


def decode_data_array(data_array, dtype):
    """
    Decode an inline base64 VTK DataArray (header_type="UInt32") back into
    a numpy array of the given dtype.
    """
    raw = base64.b64decode(data_array.text.strip())
    # '<I' = little-endian unsigned int (4 bytes), matching the UInt32
    # byte-count header VTK.write_file() prepends before the payload.
    n_bytes = struct.unpack('<I', raw[:4])[0]
    payload = raw[4:4 + n_bytes]
    return np.frombuffer(payload, dtype=dtype)


def test_save_scalar_field(tmpdir):
    mesh = CuboidMesh(4, 3, 2, 4, 3, 2)
    s = scalar_field(mesh, lambda r: r[0] + r[1] + r[2])
    vtk = VTK(mesh, directory=str(tmpdir), filename="save_scalar")
    vtk.save_scalar(s, name="s")
    path = vtk.write_file()

    assert path.endswith(".vti")
    root = ET.parse(path).getroot()
    data_array = root.find(".//DataArray[@Name='s']")
    values = decode_data_array(data_array, np.float32)
    np.testing.assert_allclose(values, s.reshape(-1).astype(np.float32))


def test_save_vector_field(tmpdir):
    mesh = CuboidMesh(4, 3, 2, 4, 3, 2)
    v = vector_field(mesh, lambda r: (r[0], r[1], r[2]))
    vtk = VTK(mesh, directory=str(tmpdir), filename="save_vector")
    vtk.save_vector(v, name="s")
    path = vtk.write_file()

    assert path.endswith(".vti")
    root = ET.parse(path).getroot()
    data_array = root.find(".//DataArray[@Name='s']")
    assert data_array.get("NumberOfComponents") == "3"
    values = decode_data_array(data_array, np.float32).reshape(-1, 3)
    np.testing.assert_allclose(values, v.reshape(-1, 3).astype(np.float32))


def test_save_scalar_field_hexagonal_mesh(tmpdir):
    mesh = HexagonalMesh(1, 3, 3)
    s = scalar_field(mesh, lambda r: r[0] + r[1])
    vtk = VTK(mesh, directory=str(tmpdir), filename="scalar_hexagonal")
    vtk.save_scalar(s, name="s")
    path = vtk.write_file()

    assert path.endswith(".vtp")
    root = ET.parse(path).getroot()

    # ".//Points/DataArray" is an ElementTree XPath: search anywhere below
    # the root (".//") for a <Points><DataArray> pair, i.e. the point
    # coordinates (the Points DataArray has no Name attribute, unlike the
    # named CellData/Polys arrays below).
    points_array = root.find(".//Points/DataArray")
    points = decode_data_array(points_array, np.float32).reshape(-1, 3)
    np.testing.assert_allclose(points, mesh.vertices.astype(np.float32))

    # [@Name='connectivity'] filters to the DataArray with that attribute,
    # since <Polys> holds two DataArrays (connectivity and offsets).
    connectivity = decode_data_array(
        root.find(".//Polys/DataArray[@Name='connectivity']"), np.int64)
    np.testing.assert_array_equal(
        connectivity.reshape(-1, 6), mesh.hexagons)

    data_array = root.find(".//CellData/DataArray[@Name='s']")
    values = decode_data_array(data_array, np.float32)
    np.testing.assert_allclose(values, s.reshape(-1).astype(np.float32))


def test_active_array_attributes_name_one_array(tmpdir):
    """
    The CellData Scalars/Vectors attributes name the array VTK makes active
    for that attribute type, so each holds a single array name. Listing every
    name space-separated matches no array and leaves nothing active, which is
    what a reader falls back on when opening the file.
    """
    mesh = CuboidMesh(nx=3, ny=2, nz=1, dx=1.0, dy=1.0, dz=1.0)
    vtk = VTK(mesh, directory=str(tmpdir), filename="active")
    vtk.save_scalar(np.zeros(mesh.n), name="Ms")
    vtk.save_scalar(np.zeros(mesh.n), name="mu_s")
    vtk.save_vector(np.zeros((mesh.n, 3)), name="spins")
    path = vtk.write_file()

    cell_data = ET.parse(path).getroot().find(".//CellData")
    assert cell_data.get("Scalars") == "Ms"
    assert cell_data.get("Vectors") == "spins"
    # the arrays that are not active are still written out
    assert {d.get("Name") for d in cell_data.findall("DataArray")} == {
        "Ms", "mu_s", "spins"}


def test_write_file_creates_nested_directory(tmpdir):
    mesh = CuboidMesh(nx=2, ny=2, nz=1, dx=1.0, dy=1.0, dz=1.0)
    vtk = VTK(mesh, directory=str(tmpdir.join("a", "b")), filename="nested")
    vtk.save_scalar(np.zeros(mesh.n), name="Ms")
    path = vtk.write_file(step=7)
    assert path.endswith(os.path.join("a", "b", "nested_000007.vti"))
    assert os.path.isfile(path)
