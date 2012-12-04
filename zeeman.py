
import numpy as np


class Zeeman(object):
    def __init__(self,H0,direction=(1,0,0)):
        self.Hx=H0*direction[0]
        self.Hy=H0*direction[1]
        self.Hz=H0*direction[2]
        
        
    def setup(self,mesh,spin):
        self.mesh=mesh
        self.dx=mesh.dx
        self.dy=mesh.dy
        self.dz=mesh.dz
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.spin=spin
        n=self.nx*self.ny*self.nz
        self.field=np.zeros((3,n))
        self.field[:,:]=[[self.Hx],[self.Hy],[self.Hz]]
        self.field.shape=(3*n)
        

    def compute_field(self):            
        return self.field
