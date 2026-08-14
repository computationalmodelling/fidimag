"""
A minimal, dependency-free writer for the modern XML-based VTK file
formats. We only ever need two dataset kinds:

    * ``.vti`` (ImageData) for the uniform CuboidMesh grid used in the
      micromagnetic and atomistic cuboid simulations.
    * ``.vtp`` (PolyData) for the HexagonalMesh, whose cells are hexagonal
      polygons defined by a list of vertices.

Field values (scalars/vectors) are always written as cell data, matching
how fidimag samples fields at cell centres. Data is stored inline as
base64-encoded binary, which keeps files compact without needing to
compute appended-data byte offsets by hand.
"""
import base64
import os
import numpy as np

from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh


def _encode(array):
    """
    Encode a numpy array as an inline base64 VTK DataArray payload: a
    leading UInt32 giving the size of the data in bytes, followed by the
    raw bytes, all base64-encoded together.
    """
    data = np.ascontiguousarray(array).tobytes()
    header = np.uint32(len(data)).tobytes()
    return base64.b64encode(header + data).decode('ascii')


def _data_array(name, array, n_components=1):
    array = np.asarray(array, dtype=np.float32)
    extra = ' NumberOfComponents="{}"'.format(n_components) if n_components > 1 else ''
    return (
        '        <DataArray type="Float32" Name="{name}"{extra} format="binary">\n'
        '          {data}\n'
        '        </DataArray>\n'
    ).format(name=name, extra=extra, data=_encode(array))


class VTK(object):
    def __init__(self, mesh, header="", directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = directory
        self.filename = filename

        if isinstance(mesh, HexagonalMesh):
            self.extension = "vtp"
        elif isinstance(mesh, CuboidMesh):
            self.extension = "vti"
        else:
            raise NotImplementedError(
                    "Mesh should be CuboidMesh or HexagonalMesh, is {}.".format(
                        mesh.__class__.__name__))

        self.header = header
        self.reset_data()

    def reset_data(self):
        self.scalars = []  # list of (name, array)
        self.vectors = []  # list of (name, array)

    def save_scalar(self, s, name="my_field", step=0):
        self.scalars.append((name, s))

    def save_vector(self, v, name="my_field", step=0):
        self.vectors.append((name, v))

    def _cell_data_block(self):
        scalar_names = " ".join(name for name, _ in self.scalars)
        vector_names = " ".join(name for name, _ in self.vectors)
        attrs = ""
        if scalar_names:
            attrs += ' Scalars="{}"'.format(scalar_names)
        if vector_names:
            attrs += ' Vectors="{}"'.format(vector_names)

        block = '      <CellData{}>\n'.format(attrs)
        for name, s in self.scalars:
            block += _data_array(name, s, n_components=1)
        for name, v in self.vectors:
            block += _data_array(name, v, n_components=3)
        block += '      </CellData>\n'
        return block

    def _write_vti(self):
        mesh = self.mesh
        nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
        extent = "0 {} 0 {} 0 {}".format(nx, ny, nz)
        origin = "{} {} {}".format(mesh.x0, mesh.y0, mesh.z0)
        spacing = "{} {} {}".format(mesh.dx, mesh.dy, mesh.dz)

        xml = '<?xml version="1.0"?>\n'
        xml += '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n'
        xml += '  <ImageData WholeExtent="{}" Origin="{}" Spacing="{}">\n'.format(
            extent, origin, spacing)
        xml += '    <Piece Extent="{}">\n'.format(extent)
        xml += self._cell_data_block()
        xml += '    </Piece>\n'
        xml += '  </ImageData>\n'
        xml += '</VTKFile>\n'
        return xml

    def _write_vtp(self):
        mesh = self.mesh
        points = np.asarray(mesh.vertices, dtype=np.float32)
        polygons = np.asarray(mesh.hexagons, dtype=np.int64)
        n_points = points.shape[0]
        n_polys = polygons.shape[0]

        connectivity = polygons.reshape(-1)
        offsets = np.arange(1, n_polys + 1, dtype=np.int64) * polygons.shape[1]

        xml = '<?xml version="1.0"?>\n'
        xml += '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n'
        xml += '  <PolyData>\n'
        xml += '    <Piece NumberOfPoints="{}" NumberOfPolys="{}">\n'.format(n_points, n_polys)
        xml += '      <Points>\n'
        xml += '        <DataArray type="Float32" NumberOfComponents="3" format="binary">\n'
        xml += '          {}\n'.format(_encode(points))
        xml += '        </DataArray>\n'
        xml += '      </Points>\n'
        xml += '      <Polys>\n'
        xml += '        <DataArray type="Int64" Name="connectivity" format="binary">\n'
        xml += '          {}\n'.format(_encode(connectivity))
        xml += '        </DataArray>\n'
        xml += '        <DataArray type="Int64" Name="offsets" format="binary">\n'
        xml += '          {}\n'.format(_encode(offsets))
        xml += '        </DataArray>\n'
        xml += '      </Polys>\n'
        xml += self._cell_data_block()
        xml += '    </Piece>\n'
        xml += '  </PolyData>\n'
        xml += '</VTKFile>\n'
        return xml

    def _render(self):
        return self._write_vti() if self.extension == "vti" else self._write_vtp()

    def write_file(self, step=0):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        filename = "{}_{:06}.{}".format(self.filename, step, self.extension)
        path = os.path.join(self.directory, filename)

        with open(path, 'w') as f:
            f.write(self._render())

        return path

    def save_as(self, path):
        """
        Write the currently stored scalar/vector data to a single file at
        `path`, without the `{filename}_{step:06d}` sequence naming used by
        `write_file`. Useful for one-off saves rather than a simulation's
        step-by-step snapshots. The mesh-appropriate extension (.vti or
        .vtp) is appended if `path` doesn't already end with it.
        """
        if not path.endswith("." + self.extension):
            path = path + "." + self.extension

        directory = os.path.dirname(path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)

        with open(path, 'w') as f:
            f.write(self._render())

        return path
