import os
import pyvtk
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
import sys


# We will try to supress the stderr output from the VTK.VtkData
# initialisations, which produce a warning when point or cell data is not
# defined at the moment of calling the class. This decorator is made for
# functions that only have variable assignments (we are not returning anything)
def ignore_stderr(func):
    def silenced(*args, **kwargs):
        with open(os.devnull, 'w') as null:
            err = sys.stderr
            sys.stderr = null
            func(*args, **kwargs)
            sys.stderr = err

    return silenced


class VTK(object):
    def __init__(self, mesh, header="", directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = directory
        self.filename = filename

        if isinstance(mesh, HexagonalMesh):
            structure = pyvtk.PolyData(points=mesh.vertices,
                                       polygons=mesh.hexagons)
        elif isinstance(mesh, CuboidMesh):
            # for keyword argument dimensions: if the mesh is made up of
            # nx * ny * nz cells, it has (nx + 1) * (ny + 1) * (nz + 1)
            # vertices.
            structure = pyvtk.RectilinearGrid(* mesh.grid)
        else:
            raise NotImplementedError(
                    "Mesh should be CuboidMesh or HexagonalMesh, is {}.".format(
                        mesh.__class__.__name__))
        self.structure = structure
        self.header = header

        self.init_VtkData(structure, header)

    @ignore_stderr
    def init_VtkData(self, structure, header):
        self.vtk_data = pyvtk.VtkData(structure, header)

    @ignore_stderr
    def reset_data(self):
        self.vtk_data = pyvtk.VtkData(self.structure, self.header)

    @ignore_stderr
    def save_scalar(self, s, name="my_field", step=0):
        self.vtk_data.cell_data.append(pyvtk.Scalars(s, name))

    def save_vector(self, v, name="my_field", step=0):
        self.vtk_data.cell_data.append(pyvtk.Vectors(v, name))

    def write_file(self, step=0):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        filename = "{}_{:06}".format(self.filename, step) + ".vtk"
        path = os.path.join(self.directory, filename)
        self.vtk_data.tofile(path)
        return path
