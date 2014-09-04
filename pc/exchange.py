import clib
from energy import Energy

class UniformExchange(Energy):
    """
    The Hamiltonian is defined as
    
        Hamiltonian = - J \sum_<i,j> S_i \cdot S_j
    
    where the brackets represent the nearest neighbours and only evaluate once
    for each pair, which means for two spins case, the total energy is -J S_1 S_2. Therefore,
    the effective field at site i is,
    
        H_i = J \sum_<i,j> S_j
    
    notice that there is no factor of 2 associated with J.
    """
    def __init__(self,J=0,name='exch'):
        self.J = J
        self.name=name
        
        self.Jx = self.J
        self.Jy = self.J
        self.Jz = self.J
        
    def compute_field(self, t=0):

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
            
        return self.field*self.mu_s_inv
            
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
