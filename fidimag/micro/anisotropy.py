import fidimag.extensions.clib
import numpy as np
from pc.constant import mu_0
from energy import Energy
import util.helper as helper

class UniaxialAnisotropy(Energy):
    """
        Ku is a number or a function
    """
    def __init__(self,Ku,axis=(1,0,0),name='anis'):
        self._Ku = Ku
        self.name = name
        self.axis = axis


    def setup(self, mesh, spin, Ms):
        super(UniaxialAnisotropy, self).setup(mesh, spin, Ms)

        self.Ms_inv_long = np.zeros(3*mesh.nxyz)

        self.Ms_inv_long.shape = (3,-1)
        for i in range(mesh.nxyz):
            if self.Ms[i] == 0.0:
                self.Ms_inv_long[:,i] = 0
            else:
                self.Ms_inv_long[:,i] = 1.0/(mu_0*self.Ms[i])

        self.Ms_inv_long.shape = (-1,)

        self.Ku = np.zeros(3*self.nxyz, dtype=np.float)
        Ku_scalar = helper.init_scalar(self._Ku, self.mesh)
        self.Ku.shape = (3,-1)
        self.Ku[0,:] = Ku_scalar[:]*self.axis[0]
        self.Ku[1,:] = Ku_scalar[:]*self.axis[1]
        self.Ku[2,:] = Ku_scalar[:]*self.axis[2]
        self.Ku.shape = (-1,)

    def compute_field(self, t=0):
        clib.compute_anisotropy(self.spin,
                                self.field,
                                self.Ku,
                                self.nxyz)

        return self.field*self.Ms_inv_long

    def compute_energy(self):
        self.energy=clib.compute_anisotropy_energy(self.spin,
                                            self.Ku,
                                            self.nxyz)

        return self.energy*self.mesh.cellsize
