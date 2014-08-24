import clib 
import numpy as np


class Demag(object):
    def __init__(self, name='demag', oommf = True):
        self.name = name
        self.oommf = oommf
        
    def setup(self, mesh, spin, Ms):
        self.mesh = mesh
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.n = mesh.nxyz
        self.field = np.zeros(3 * self.n, dtype=np.float)
        
        self.mu_s = np.zeros(mesh.nxyz,dtype=np.float)
        
        self.mu_s[:] = Ms

        self.demag=clib.FFTDemag(self.dx,self.dy,self.dz,
                                 self.nx,self.ny,self.nz,
                                 self.oommf)
        
    def compute_field(self, t=0):
        self.demag.compute_field(self.spin,self.mu_s,self.field)
        return self.field

    def compute_energy(self):
        
        #energy=self.demag.compute_energy(self.spin,self.mu_s,self.field)
        
        return 0.0
    
