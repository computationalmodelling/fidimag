import os
import numpy as np
from fidimag.common.vtk import VTK


class SaveVTK():

    def __init__(self, mesh, name='unnamed'):
        self.name = name

        # Initiate a VTK object
        self.VTK = VTK(mesh,
                       directory='{}_vtks'.format(name),
                       filename=name
                       )

    def save_vtk(self, m1, Ms, step=0, vtkname='m'):

        # Here we save both Ms and spins as cell data
        self.VTK.save_scalar(Ms, name='Ms')
        self.VTK.save_vector(m1, name='spins')

        # Now we set the name of the file to start with m_ as default
        self.VTK.filename = vtkname
        self.VTK.write_file(step=step)

    def save_vtk_scalar(self, scalar_data_array, step=0, vtkname='skx'):
        # This saves the skyrmion number or any scalar data
        # CHECK: that the array skx_num is being passed with the correct shape
        # TODO: put more generic names as default names
        self.VTK.save_scalar(scalar_data_array, name='skx_num')

        # These files will be saved starting with: skx_ or other name
        # and saved to a different folder instead of {}_vtks
        self.VTK.filename = vtkname
        self.VTK.directory = self.name + '_' + vtkname
        self.VTK.write_file(step=step)
