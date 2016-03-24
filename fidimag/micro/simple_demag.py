import numpy as np


class SimpleDemag(object):
    """
    Demagnetising field for thin films in the i-direction.
    Hj = Hk = 0 and Hi = - Mi.

    """
    def __init__(self, Nx=0, Ny=0.5, Nz=0.5, name='SimpleDemag'):
        """
        field_strength is Ms by default

        """
        
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        
        self.name = name

    def setup(self, mesh, spin, Ms):
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n

        self.Ms = Ms
        self.Ms_long = np.zeros(3 * mesh.n)

        self.Ms_long.shape = (3, -1)
        for i in range(mesh.n):
            self.Ms_long[:, i] = Ms[i]

        self.Ms_long.shape = (-1,)
        self.field = np.zeros(3 * self.n)
        #self.field[:] = helper.init_vector(self.H0, self.mesh)


    def compute_field(self, t=0):
        self.spin.shape = (3,-1)
        self.field.shape = (3,-1)
        self.field[0][:] = -self.Nx*self.spin[0][:]*self.Ms[:]
        self.field[1][:] = -self.Ny*self.spin[1][:]*self.Ms[:]
        self.field[2][:] = -self.Nz*self.spin[2][:]*self.Ms[:]
        self.spin.shape = (-1,)
        self.field.shape = (-1,)
        return self.field

    def compute_energy(self):

        mu_0 = 4*np.pi*1e-7

        sf = 0.5* self.field * self.spin * self.Ms_long * mu_0

        energy = -np.sum(sf)

        return energy * self.mesh.cellsize
