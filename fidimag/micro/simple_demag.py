import numpy as np
from .energy import Energy


class SimpleDemag(Energy):
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
        sf = -0.5 * self.field * self.spin * mu_0
        energy = np.sum(sf.reshape(-1, 3), axis=1) * self.Ms

        return energy * self.mesh.cellsize
