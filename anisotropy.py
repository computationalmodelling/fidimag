import clib 
import numpy as np


class Anisotropy(object):
    def __init__(self,K,direction=(1,0,0)):
        self.Kx=K*direction[0]
        self.Ky=K*direction[1]
        self.Kz=K*direction[2]
        
        
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
        self.field=np.zeros(3*n)

    def compute_field(self):
        clib.compute_anisotropy(self.spin,
                                self.field,
                                self.Kx,
                                self.Ky,
                                self.Kz,
                                self.nx,
                                self.ny,
                                self.nz)
                                      
        return self.field
