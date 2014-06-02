import clib 
import numpy as np
import pccp.util.helper as helper

class Anisotropy(object):
    """
    Hamiltonian = -Ku_x*S_x^2 - Ku_y*S_y^2 -Ku_z*S_z^2 
    ==>H=2D*S_x
    """
    def __init__(self,Ku,name='anis'):
        self._Ku = Ku
        self.name=name
    
    
    def setup(self,mesh,spin,mu_s_inv=1.0, pbc=None):
        self.mesh=mesh
        self.spin=spin
            
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)
        self.energy=0
        self.mu_s_inv = mu_s_inv
        
        self.Ku = np.zeros(3*self.nxyz)
        self.Ku[:] = helper.init_vector(self._Ku, self.mesh)

    def compute_field(self, t=0):
        clib.compute_anisotropy(self.spin,
                                self.field,
                                self.Ku,
                                self.nxyz)
                                      
        return self.field*self.mu_s_inv
    
    def compute_energy(self):
        self.energy=clib.compute_anisotropy_energy(self.spin,
                                            self.Ku,
                                            self.nxyz)
        
        return self.energy
