import os
import pyvtk
import numpy as np


class SaveVTK():

    def __init__(self, mesh, name='unnamed'):
        self.mesh = mesh
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.name = '%s_vtks' % name
        xyz = mesh.coordinates
        self.x = np.array(xyz[:, 0], dtype='float32')
        self.y = np.array(xyz[:, 1], dtype='float32')
        self.z = np.array(xyz[:, 2], dtype='float32')

        # build a new index since we have used difference order
        ids = [self.mesh.index(i, j, k)
               for i in range(self.nx)
               for j in range(self.ny)
               for k in range(self.nz)]

        self.ids = np.array(ids)

        self.pos = [self.mesh.coordinates[self.ids[i]] for i in range(mesh.n)]

    def save_vtk(self, m1, step=0, vtkname='m'):

        if not os.path.exists(self.name):
            os.makedirs(self.name)

        pos = pyvtk.StructuredGrid([self.nx, self.ny, self.nz], self.pos)
        m = m1.copy()
        m.shape = (-1, 3)
        data = pyvtk.PointData(pyvtk.Vectors(m[self.ids], 'm'))
        m.shape = (-1,)

        vtk = pyvtk.VtkData(pos, data, 'spins')

        vtk.tofile("%s/%s.%06d" % (self.name, vtkname, step), 'binary')

    def save_vtk_scalar(self, skx_num, step=0, vtkname='skx'):

        folder_name = self.name + '_' + vtkname
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        pos = pyvtk.StructuredGrid([self.nx, self.ny, self.nz], self.pos)

        data = pyvtk.PointData(
            pyvtk.Scalars(skx_num[self.ids], 'skx_num', lookup_table='default'))

        vtk = pyvtk.VtkData(pos, data, 'skx_num')

        vtk.tofile("%s/%s.%06d" % (folder_name, vtkname, step), 'binary')


class SaveVTK_unstructured():

    def __init__(self, sim, name='unnamed'):
        """

        Saves VTK files with non-zero Ms points

        The current problem of the UnstructuredGrid is that the surface
        cannot be visualised in Mayavi or Paraview, but
        it saves some memory when not considering points with Ms = 0

        Needs:

            simulation --> to get the mesh and
                           to mask the values where Ms is zero

        """

        self.mesh = sim.mesh
        # Number of elements in the corresponding directions
        self.nx = self.mesh.nx
        self.ny = self.mesh.ny
        self.nz = self.mesh.nz
        # Corresponding finite differences sizes
        self.dx = self.mesh.dx
        self.dy = self.mesh.dy
        self.dz = self.mesh.dz

        self.name = '%s_vtks' % name

        xyz = np.array(self.mesh.coordinates)

        self.x = np.array(xyz[:, 0], dtype='float32')
        self.y = np.array(xyz[:, 1], dtype='float32')
        self.z = np.array(xyz[:, 2], dtype='float32')

        # build a new index since we have used difference order
        ids = [self.mesh.index(i, j, k)
               for k in range(self.nz)
               for j in range(self.ny)
               for i in range(self.nx)]

        self.ids = np.array(ids)

        self.pos = []
        for i in range(len(ids)):
            self.pos.append(self.mesh.coordinates[self.ids[i]])

        # Reorder the saturation magnetisation values
        self.mask = sim.Ms[self.ids]
        # This will return a masked array and a mask array
        # with True where the values are zero
        self.mask = np.ma.masked_values(self.mask, 0)
        # Here we call the mask from the masked object
        # and invert the mask to get only the nonzero Ms values
        self.mask = ~self.mask.mask

    def save_vtk(self, m, step=0, vtkname='m'):

        if not os.path.exists(self.name):
            os.makedirs(self.name)

        pos = pyvtk.UnstructuredGrid(np.array(self.pos)[self.mask])

        m.shape = (3, -1)

        data = pyvtk.PointData(pyvtk.Vectors(np.transpose(m)[self.ids][self.mask],
                                             'm')
                               )
        m.shape = (-1, )

        vtk = pyvtk.VtkData(pos, data, 'spins')

        vtk.tofile("%s/%s.%06d" % (self.name, vtkname, step), 'binary')

    def save_vtk_scalar(self, skx_num, step=0, vtkname='skx'):

        folder_name = self.name + '_' + vtkname
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        pos = pyvtk.UnstructuredGrid(self.pos[self.mask])

        data = pyvtk.PointData(pyvtk.Scalars(skx_num[self.ids][self.mask],
                                             'skx_num',
                                             lookup_table='default'))

        vtk = pyvtk.VtkData(pos, data, 'skx_num')

        vtk.tofile("%s/%s.%06d" % (folder_name, vtkname, step), 'binary')
