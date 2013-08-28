import clib 
import numpy as np


class UniformExchange(object):
    """
    The Hamiltonian is defined as
    
        Hamiltonian = -J \sum_<i,j> S_i \cdot S_j
    
    where the brackets represent the nearest neighbours and only evaluate once
    for each pair, which means for two spins case, the total energy is -J S_1 S_2. Therefore,
    the effective field at site i is,
    
        H_i = J \sum_<i,j> S_j
    
    notice that there is no factor of 2 associated with J.
    """
    def __init__(self,J,name='exch'):
        self.J=J
        self.name=name
        
    def setup(self,mesh,spin,unit_length=1.0,mu_s=1.0):
        self.J_mu_s=self.J/mu_s
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
                                      self.J_mu_s,
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