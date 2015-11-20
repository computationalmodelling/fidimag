import os
import numpy as np
from fidimag.common.vtk import VTK


class SaveVTK():

    def __init__(self, mesh, name='unnamed'):
        self.name = name
        self.VTK = VTK(mesh,
                       directory='{}_vtks'.format(name),
                       filename=name
                       )

    def save_vtk(self, m1, Ms, step=0, vtkname='m'):

        self.VTK.save_scalar(Ms, name='Ms')
        self.VTK.save_vector(m1, name='spins')

        self.VTK.filename = vtkname
        self.VTK.write_file(step=step)

    def save_vtk_scalar(self, skx_num, step=0, vtkname='skx'):

        self.VTK.save_scalar(skx_num, name='skx_num')

        self.VTK.filename = vtkname
        self.VTK.directory = self.name + '_' + vtkname
        self.VTK.write_file(step=step)
