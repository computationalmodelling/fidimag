import clib 
import numpy as np
from constant import mu_0

class UniformExchange(object):
    """
    The Hamiltonian is defined as
    
        Hamiltonian = - J \sum_<i,j> S_i \cdot S_j
    
    where the brackets represent the nearest neighbours and only evaluate once
    for each pair, which means for two spins case, the total energy is -J S_1 S_2. Therefore,
    the effective field at site i is,
    
        H_i = J \sum_<i,j> S_j
    
    notice that there is no factor of 2 associated with J.
    """
    def __init__(self,J=0,name='exch', A=None):
        self.J = J
        self.A = A
        self.name=name
        
        self.Jx = self.J
        self.Jy = self.J
        self.Jz = self.J
        
    def setup(self,mesh,spin,mu_s_inv,pbc=None):
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
        self.energy=np.zeros(3*n)
        self.total_energy=0
        self.pbc = pbc
        self.mu_s_inv = mu_s_inv
                
        if self.A is not None:
            self.mu_s_inv[:] *= (2.0/mu_0)/(mesh.unit_length**2)
            print 'step1',self.mu_s_inv
    
        self.xperiodic = 0
        self.yperiodic = 0
        
        if pbc=='1d':
            self.xperiodic = 1
        elif pbc=='2d':
            self.xperiodic = 1
            self.yperiodic = 1

    def compute_field(self, t=0):
        if self.A is not None:
            clib.compute_exchange_field_c(self.spin,
                                      self.field,
                                      self.energy,
                                      self.A,
                                      self.dx,
                                      self.dy,
                                      self.dz,
                                      self.nx,
                                      self.ny,
                                      self.nz,
                                      self.xperiodic,
                                      self.yperiodic)
            
        else:
            clib.compute_exchange_field(self.spin,
                                      self.field,
                                      self.energy,
                                      self.Jx,
                                      self.Jy,
                                      self.Jz,
                                      self.nx,
                                      self.ny,
                                      self.nz,
                                      self.xperiodic,
                                      self.yperiodic)
        
        self.total_energy = np.sum(self.energy)/2.0
        
        print 'hahahah',self.mu_s_inv
        return self.field*self.mu_s_inv
        
    def compute_energy(self):
        
        #since we are not always calling this function, so it's okay to call compute_field again
        self.compute_field()
        
        self.total_energy = np.sum(self.energy)/2.0
        
        return self.total_energy
    
    def compute_energy_directly(self):
        
        energy=clib.compute_exchange_energy(self.spin,
                                                 self.Jx,
                                                 self.Jy,
                                                 self.Jz,
                                                 self.nx,    
                                                 self.ny,
                                                 self.nz,
                                                 self.xperiodic,
                                                 self.yperiodic)
                                                 
        return energy
