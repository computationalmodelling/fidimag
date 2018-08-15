import fidimag.extensions.dipolar as clib
import numpy as np
import os
from .energy import Energy

mu_0 = 4 * np.pi * 1e-7


default_options={
    'pbc_2d_error': 1e-10,
    'sample_repeat_nx': -1,
    'sample_repeat_ny': -1,
    'tensor_file_name': '2dpbc_tensors',
}


class Demag(Energy):

    def __init__(self, name='Demag', pbc_2d=False,
                 pbc_options=default_options):
        self.name = name
        self.oommf = True
        self.pbc_2d = pbc_2d
        self.pbc_options = pbc_options
        self.jac = False

    def setup(self, mesh, spin, Ms, Ms_inv):
        super(Demag, self).setup(mesh, spin, Ms, Ms_inv)

        if self.pbc_2d is True:

            self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                       self.nx, self.ny, self.nz, tensor_type='2d_pbc')

            nxyz = self.nx*self.ny*self.nz
            tensors = np.zeros(6*nxyz, dtype=np.float)

            pbc_2d_error = 1e-10
            sample_repeat_nx = -1
            sample_repeat_ny = -1
            asymptotic_radius = 32.0
            dipolar_radius = 10000.0
            tensor_file_name = ''

            options = self.pbc_options

            if 'sample_repeat_nx' in options:
                sample_repeat_nx = options['sample_repeat_nx']
                sample_repeat_ny = options['sample_repeat_ny']
            if 'relative_tensor_error' in options:
                pbc_2d_error = options['relative_tensor_error']
            if 'asymptotic_radius' in options:
                asymptotic_radius = options['asymptotic_radius']
            if 'dipolar_radius' in options:
                dipolar_radius = options['dipolar_radius']
            if 'tensor_file_name' in options:
                tensor_file_name = options['tensor_file_name']

            if len(tensor_file_name) > 0 :
                if not (os.path.exists(tensor_file_name+'.npz')):

                    self.demag.compute_tensors_2dpbc(tensors, pbc_2d_error, sample_repeat_nx, sample_repeat_ny, dipolar_radius)
                else:
                    npzfile = np.load(tensor_file_name+'.npz')
                    geo_info = npzfile['geo']
                    geo_arr = np.array([self.nx,self.ny, self.nz, self.dx, self.dy, self.dz], dtype=np.float)
                    #print 'info', geo_info
                    if not np.allclose(geo_arr, geo_info):
                        self.demag.compute_tensors_2dpbc(tensors, pbc_2d_error, sample_repeat_nx, sample_repeat_ny, dipolar_radius)
                    else:
                        tensors = npzfile['tensors']

                geo_arr = np.array([self.nx, self.ny, self.nz, self.dx, self.dy, self.dz], dtype=np.float)
                np.savez(tensor_file_name+'.npz', geo=geo_arr, tensors=tensors)

            else:
                self.demag.compute_tensors_2dpbc(tensors, pbc_2d_error, sample_repeat_nx, sample_repeat_ny, dipolar_radius)

            #print tensors
            self.demag.fill_demag_tensors(tensors)


        else:
            self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                       self.nx, self.ny, self.nz,
                                       tensor_type='demag')



    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        self.demag.compute_field(m, self.Ms, self.field)

        return self.field

    def compute_exact(self):
        field = np.zeros(3 * self.mesh.n)
        self.demag.compute_exact(self.spin, self.Ms, field)
        return field

    def compute_energy(self):

        self.compute_field()
        energy = self.demag.compute_energy(self.spin, self.Ms,
                                           self.field, self.energy)

        self.energy *= mu_0 * (self.mesh.dx *
                               self.mesh.dy *
                               self.mesh.dz *
                               self.mesh.unit_length ** 3.)

        return energy * mu_0 * (self.mesh.dx *
                                self.mesh.dy *
                                self.mesh.dz *
                                self.mesh.unit_length ** 3.)
