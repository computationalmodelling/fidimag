import fidimag.extensions.dipolar as clib
import numpy as np

mu_0 = 4 * np.pi * 1e-7


default_options={
    'pbc_2d_error':1e-10, 
    'sample_repeat_nx':-1,
    'sample_repeat_ny':-1,
}

class Demag(object):

    def __init__(self, name='demag', pbc_2d=False, pbc_options=default_options):
        self.name = name
        self.oommf = True
        self.pbc_2d = pbc_2d
        self.pbc_options = pbc_options

    def setup(self, mesh, spin, Ms):
        self.mesh = mesh
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.field = np.zeros(3 * mesh.n, dtype=np.float)

        self.Ms = Ms

        if self.pbc_2d is True:
            self.demag = clib.FFTDemag(self.dx, self.dy, self.dz, 
                                       self.nx, self.ny, self.nz, tensor_type='2d_pbc', options=self.pbc_options)

        else:
            self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                       self.nx, self.ny, self.nz, tensor_type='demag')



    def compute_field(self, t=0):

        self.demag.compute_field(self.spin, self.Ms, self.field)

        return self.field

    def compute_exact(self):
        field = np.zeros(3 * self.mesh.n)
        self.demag.compute_exact(self.spin, self.Ms, field)
        return field

    def compute_energy(self):

        self.compute_field()
        energy = self.demag.compute_energy(self.spin, self.Ms, self.field)

        return energy * mu_0 * (self.mesh.dx *
                                self.mesh.dy *
                                self.mesh.dz *
                                self.mesh.unit_length ** 3.)
