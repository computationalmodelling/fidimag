import os
from tvtk.api import tvtk, write_data  # distributed with enthought/mayavi2
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh


class VTK(object):
    def __init__(self, mesh, directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = directory
        self.filename = filename

        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        if isinstance(mesh, HexagonalMesh):
            self.grid = tvtk.PolyData(points=mesh.vertices, polys=mesh.hexagons)
        elif isinstance(mesh, CuboidMesh):
            # if the mesh is made up of nx * ny * nz cells, it has
            # (nx + 1) * (ny + 1) * (nz + 1) vertices.
            self.grid = tvtk.StructuredGrid(dimensions=[ni + 1 for ni in mesh.size])
            self.grid.points = mesh.grid
        else:
            raise NotImplementedError("Mesh should be CuboidMesh or HexagonalMesh.")

    def save_scalar(self, s, step=0, name="my_field"):
        self.grid.cell_data.scalars = s.copy()
        self.grid.cell_data.scalars.name = name

        filename = "{}_{:06}".format(self.filename, step)
        write_data(self.grid, os.path.join(self.directory, filename))

    def save_vector(self, v, step=0, name="my_field"):
        self.grid.cell_data.vectors = v.copy()
        self.grid.cell_data.vectors.name = name

        filename = "{}_{:06}".format(self.filename, step)
        write_data(self.grid, os.path.join(self.directory, filename))
