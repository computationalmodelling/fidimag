import os
from tvtk.api import tvtk, write_data  # distributed with enthought/mayavi2
# if the import above fails due to missing X display, export ETS_TOOLKIT='null'
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh


class VTK(object):
    def __init__(self, mesh, directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = directory
        self.filename = filename

        if isinstance(mesh, HexagonalMesh):
            self.grid = tvtk.PolyData(points=mesh.vertices, polys=mesh.hexagons)
        elif isinstance(mesh, CuboidMesh):
            # if the mesh is made up of nx * ny * nz cells, it has
            # (nx + 1) * (ny + 1) * (nz + 1) vertices.
            self.grid = tvtk.StructuredGrid(dimensions=[ni + 1 for ni in mesh.size])
            self.grid.points = mesh.grid
        else:
            raise NotImplementedError("Mesh should be CuboidMesh or HexagonalMesh.")

    def save_scalar(self, s, name="my_field", step=0):
        self.grid.cell_data.scalars = s.copy()
        self.grid.cell_data.scalars.name = name

    def save_vector(self, v, name="my_field", step=0):
        self.grid.cell_data.vectors = v.copy()
        self.grid.cell_data.vectors.name = name

    def write_file(self, step=0):
        # We will create the directory only when writing the file
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        filename = "{}_{:06}".format(self.filename, step)
        write_data(self.grid, os.path.join(self.directory, filename))
