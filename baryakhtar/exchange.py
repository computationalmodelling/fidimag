import baryakhtar_clib as clib 
import numpy as np

class UniformExchange(object):
    """
    compute the exchange field.
    """
    def __init__(self, A, chi, name='exch'):
        self.A = A
        self.chi = chi
        self.name = name    
    
    def setup(self,mesh,m,Me):
        self.mesh=mesh
        self.dx=mesh.dx*mesh.unit_length
        self.dy=mesh.dy*mesh.unit_length
        self.dz=mesh.dz*mesh.unit_length
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.spin = m
        self.Me = Me
        self.field = np.zeros(3*mesh.nxyz)
        
        if self.chi == 0.0:
            self.chi_inv = 0
        else:
            self.chi_inv = 1.0/self.chi


    def compute_field(self, t=0):
        
        clib.compute_exchange_field_baryakhtar(self.spin,
                                      self.field,
                                      self.Me,
                                      self.chi_inv,
                                      self.A,
                                      self.dx,
                                      self.dy,
                                      self.dz,
                                      self.nx,
                                      self.ny,
                                      self.nz)
        
        return self.field
        
    def compute_energy(self):
        
        return 0.0
