import clib 
import numpy as np


class Anisotropy(object):
    """
    Hamiltonian=-D*S_x^2 ==>H=2D*S_x
    """
    def __init__(self,D,direction=(1,0,0),name='anis'):
        self.Dx=D*direction[0]
        self.Dy=D*direction[1]
        self.Dz=D*direction[2]
        self.name=name
    
    
    def setup(self,mesh,spin,mu_s_inv=1.0, pbc=None):
        self.mesh=mesh
        self.spin=spin
        
        #TODO: change Dx, Dy, Dz to arrays
        self.Dx=self.Dx
        self.Dy=self.Dy
        self.Dz=self.Dz
        
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)
        self.energy=0
        self.mu_s_inv = mu_s_inv

    def compute_field(self, t=0):
        clib.compute_anisotropy(self.spin,
                                self.field,
                                self.Dx_mu_s,
                                self.Dy_mu_s,
                                self.Dz_mu_s,
                                self.nxyz)
                                      
        return self.field*self.mu_s_inv
    
    def compute_energy(self):
        self.energy=clib.compute_anisotropy_energy(self.spin,
                                       self.Dx,
                                       self.Dy,
                                       self.Dz,
                                       self.nxyz)
        
        return self.energy
