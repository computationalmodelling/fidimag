import clib 
import numpy as np


class Anisotropy(object):
    """
    Hamiltonian=-D*S_x^2 ==>H=2D*S_x
    """
    def __init__(self,D,direction=(1,0,0),mu_s=1.0):
        
        self.Dx=D*direction[0]/mu_s
        self.Dy=D*direction[1]/mu_s
        self.Dz=D*direction[2]/mu_s
        
        
    def setup(self,mesh,spin,unit_length=1.0):
        self.mesh=mesh
        self.spin=spin
        self.nxyz=mesh.nxyz
        self.field=np.zeros(3*self.nxyz)

    def compute_field(self):
        clib.compute_anisotropy(self.spin,
                                self.field,
                                self.Dx,
                                self.Dy,
                                self.Dz,
                                self.nxyz)
                                      
        return self.field
