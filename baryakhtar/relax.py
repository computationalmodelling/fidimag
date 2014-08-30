import baryakhtar_clib as clib 
import numpy as np

class Relaxation(object):
    """
    compute the relaxation field related to exchange field.
    """
    def __init__(self, chi, name='relax'):
        self.chi = chi
        self.name = name    
    
    def setup(self, mesh, m, Ms):
        self.mesh=mesh
        self.dx=mesh.dx*mesh.unit_length
        self.dy=mesh.dy*mesh.unit_length
        self.dz=mesh.dz*mesh.unit_length
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.spin = m
        self.Ms = Ms
        self.field = np.zeros(3*mesh.nxyz)
        
        if self.chi == 0.0:
            self.chi_inv = 0
        else:
            self.chi_inv = 1.0/self.chi


    def compute_field(self, t=0):
        
        clib.compute_relaxation_field(self.spin,
                                      self.field,
                                      self.Ms,
                                      self.chi_inv,
                                      self.mesh.nxyz)
        
        return self.field
        
    def compute_energy(self):
        
        return 0.0


class Laplace(object):
    """
        compute the laplace for given field.
    """
    def __init__(self, mesh):
        self.dx=mesh.dx*mesh.unit_length
        self.dy=mesh.dy*mesh.unit_length
        self.dz=mesh.dz*mesh.unit_length
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.field = np.zeros(3*mesh.nxyz)
    
    
    def compute_laplace_field(self, h):
        
        clib.compute_laplace_field(h, self.field,
                                   self.dx,
                                   self.dy,
                                   self.dz,
                                   self.nx,
                                   self.ny,
                                   self.nz)
                                   
        return self.field
