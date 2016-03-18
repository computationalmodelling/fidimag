import os
import pyvtk 
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh


class VTK(object):
    def __init__(self, mesh, header="", directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = directory
        self.filename = filename

        if isinstance(mesh, HexagonalMesh):
            structure = pyvtk.PolyData(points=mesh.vertices, polygons=mesh.hexagons)
            self.vtk_data = pyvtk.VtkData(structure, header)
        elif isinstance(mesh, CuboidMesh):
            # for keyword argument dimensions: if the mesh is made up of
            # nx * ny * nz cells, it has (nx + 1) * (ny + 1) * (nz + 1)
            # vertices.
            structure = pyvtk.RectilinearGrid(* mesh.grid)
            self.vtk_data = pyvtk.VtkData(structure, header)
        else:
            raise NotImplementedError(
                    "Mesh should be CuboidMesh or HexagonalMesh, is {}.".format(
                        mesh.__class__.__name__))

    def save_scalar(self, s, name="my_field", step=0):
        self.vtk_data.cell_data.append(pyvtk.Scalars(s, name))

    def save_vector(self, v, name="my_field", step=0):
        self.vtk_data.cell_data.append(pyvtk.Vectors(v, name))

    def write_file(self, step=0):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        filename = "{}_{:06}".format(self.filename, step)
        path = os.path.join(self.directory, filename)
        self.vtk_data.tofile(path)
