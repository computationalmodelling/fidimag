import os
import pyvtk
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
import sys

log = pyvtk.logging.getLogger(pyvtk.__name__)
log.setLevel(pyvtk.logging.ERROR)


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

    def init_VtkData(self, structure, header):
        self.vtk_data = pyvtk.VtkData(structure, header)

    def reset_data(self):
        self.vtk_data = pyvtk.VtkData(self.structure, self.header)

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

# Attempt to get rid of old pyvtk
import vtk
import vtk.util.numpy_support as nps
from pathlib import Path

class VTKImageData(object):
    def __init__(self, mesh, header="", directory=".", filename="unnamed"):
        self.mesh = mesh
        self.directory = Path(directory)
        self.vtkArrays = []

        if isinstance(mesh, HexagonalMesh):
            # structure = pyvtk.PolyData(points=mesh.vertices,
            #                            polygons=mesh.hexagons)
            # structure = 
            raise NotImplementedError('')
        elif isinstance(mesh, CuboidMesh):
            self.filename = Path(filename).with_suffix('.vti')
            # for keyword argument dimensions: if the mesh is made up of
            # nx * ny * nz cells, it has (nx + 1) * (ny + 1) * (nz + 1)
            # vertices.
            # structure = pyvtk.RectilinearGrid(* mesh.grid)
            structure = vtk.vtkImageData()
            structure.SetDimensions(mesh.nx + 1, mesh.ny + 1, mesh.nz + 1)
            structure.SetSpacing(mesh.dx, mesh.dy, mesh.dz)
            structure.SetOrigin(mesh.x0, mesh.y0, mesh.z0)

            self.writer = vtk.vtkXMLImageDataWriter()
            self.writer.SetFileName(filename)
            self.writer.SetInputData(structure)
        else:
            raise NotImplementedError(
                    "Mesh should be CuboidMesh or HexagonalMesh, is {}.".format(
                        mesh.__class__.__name__))

        self.structure = structure
        self.structureCellData = self.structure.GetCellData()
        self.header = header

    def reset_data(self):
        for idx in self.vtkArrays:
            self.structureCellData.RemoveArray(idx)
            self.structure.Modified()
        # self.vtk_data = pyvtk.VtkData(self.structure, self.header)

    def save_array(self, a, name="my_field"):
        """Vectors must be reshaped, e.g. (mesh.n, 3)"""
        # self.vtk_data.cell_data.append(pyvtk.Scalars(s, name))
        # TODO: add type, double or float to vtk array
        # Deep copy, check: https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        mData = nps.numpy_to_vtk(a, deep=True)
        mData.SetName(name)
        self.structureCellData.AddArray(mData)
        self.structure.Modified()

    def write_file(self, step=0):
        self.directory.mkdir(exist_ok=True)
        path = self.directory / f"{self.filename.stem}_{step:06}.vti"
        self.writer.SetFileName(str(path))
        self.writer.Write()
        return path
