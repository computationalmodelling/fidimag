import clib 
import numpy as np


class Demag(object):
    def __init__(self, name='demag', oommf=False):
        self.name = name
        self.oommf = oommf
        
    def setup(self, mesh, spin, mu_s_inv=1, pbc=None):
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
        unit_length = mesh.unit_length
        self.mu_s = np.zeros(mesh.nxyz,dtype=np.float)
        
        #note that the 1e-7 comes from \frac{\mu_0}{4\pi}
        self.scale = 1e-7 / unit_length**3
        
        for i in range(self.n):
            if mu_s_inv[i] == 0.0:
                self.mu_s[i] = 0.0
            else: 
                if self.oommf:
                    self.mu_s[i] = 1.0/mu_s_inv[i]
                else:
                    self.mu_s[i] = 1.0/mu_s_inv[i]*self.scale
        
        
        self.demag=clib.FFTDemag(self.dx,self.dy,self.dz,
                                 self.nx,self.ny,self.nz,
                                 self.oommf)
        
    def compute_field(self, t=0):
        self.demag.compute_field(self.spin,self.mu_s,self.field)
        return self.field
    
    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin,self.mu_s,field)
        return field

    def compute_energy(self):
        
        energy=self.demag.compute_energy(self.spin,self.mu_s,self.field)
        
        return energy/self.scale
    
