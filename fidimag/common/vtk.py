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
    data = array if isinstance(array, bytes) else np.ascontiguousarray(array).tobytes()
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


def _string_array(name, text):
    """
    A VTK String DataArray holding one value. VTK stores strings as a run of
    null-terminated bytes, so the payload is the text plus a terminator,
    framed and base64-encoded exactly like a numeric array. Going through
    base64 also means the text needs no XML escaping.
    """
    return (
        '      <DataArray type="String" Name="{name}" NumberOfTuples="1" format="binary">\n'
        '        {data}\n'
        '      </DataArray>\n'
    ).format(name=name, data=_encode(text.encode('utf-8') + b'\x00'))


def _fidimag_version():
    try:
        from importlib.metadata import version
        return version('fidimag')
    except Exception:
        return 'unknown'


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

    def _default_header(self):
        """
        A one-line description of what wrote the file and of the mesh it was
        written on, so a .vti/.vtp is still self-describing once it is sitting
        in a vtks/ directory months later.
        """
        mesh = self.mesh
        if isinstance(mesh, HexagonalMesh):
            # A HexagonalMesh is two-dimensional. Its nz and dz exist only as
            # placeholders (both 1), so report what actually defines it: the
            # two lattice counts, the hexagon radius and the alignment.
            geometry = "HexagonalMesh nx={} ny={} n={} radius={:g} alignment={}".format(
                mesh.nx, mesh.ny, mesh.n, mesh.radius, mesh.alignment)
        else:
            geometry = "CuboidMesh nx={} ny={} nz={} n={} dx={:g} dy={:g} dz={:g}".format(
                mesh.nx, mesh.ny, mesh.nz, mesh.n, mesh.dx, mesh.dy, mesh.dz)

        return "fidimag {} | {} | unit_length={:g}".format(
            _fidimag_version(), geometry, mesh.unit_length)

    def _field_data_block(self):
        """
        Dataset-level metadata, which readers show alongside the file rather
        than as a field on the mesh. The legacy .vtk format had a title line
        for this; the XML formats use FieldData instead.
        """
        block = '    <FieldData>\n'
        block += _string_array('fidimag', self._default_header())
        if self.header:
            block += _string_array('description', self.header)
        block += '    </FieldData>\n'
        return block

    def reset_data(self):
        self.scalars = []  # list of (name, array)
        self.vectors = []  # list of (name, array)

    def save_scalar(self, s, name="my_field", step=0):
        self.scalars.append((name, s))

    def save_vector(self, v, name="my_field", step=0):
        self.vectors.append((name, v))

    def _cell_data_block(self):
        # These attributes name the array VTK should make active for each
        # attribute type, so each takes a single name. Listing several
        # space-separated matches no array at all and leaves nothing active,
        # so name the first of each and let the rest simply be present.
        attrs = ""
        if self.scalars:
            attrs += ' Scalars="{}"'.format(self.scalars[0][0])
        if self.vectors:
            attrs += ' Vectors="{}"'.format(self.vectors[0][0])

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
        xml += self._field_data_block()
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
        xml += self._field_data_block()
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

    def _write(self, path):
        directory = os.path.dirname(path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)

        with open(path, 'w') as f:
            f.write(self._render())

        return path

    def write_file(self, step=0):
        filename = "{}_{:06}.{}".format(self.filename, step, self.extension)
        return self._write(os.path.join(self.directory, filename))

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

        return self._write(path)
