import clib 
import numpy as np


class UniformExchange(object):
    def __init__(self,J,name='exch'):
        self.J=J
        self.name=name
        
    def setup(self,mesh,spin,unit_length=1.0,mu_s=1.0):
        self.J/=mu_s
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
        self.energy=0

    def compute_field(self):
        clib.compute_uniform_exchange(self.spin,
                                      self.field,
                                      self.J,
                                      self.dx,
                                      self.dy,
                                      self.dz,
                                      self.nx,
                                      self.ny,
                                      self.nz)                     
        return self.field

    def compute_energy(self):
        
        self.energy=clib.compute_exchange_energy(self.spin,
                                                 self.J,
                                                 self.nx,    
                                                 self.ny,
                                                 self.nz)
        return self.energy