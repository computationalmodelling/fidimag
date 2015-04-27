import fidimag.extensions.micro_clib as micro_clib
from energy import Energy
#from constant import mu_0

class UniformExchange(Energy):
    """
    Compute the exchange field in micromagnetics.
    """
    def __init__(self, A, name='exch'):
        self.A = A
        self.name = name
        
    def compute_field(self, t=0):

        micro_clib.compute_exchange_field_micro(self.spin,
                                                self.field,
                                                self.energy,
                                                self.Ms_inv,
                                                self.A,
                                                self.dx,
                                                self.dy,
                                                self.dz,
                                                self.nx,
                                                self.ny,
                                                self.nz,
                                                self.xperiodic,
                                                self.yperiodic)
        
        
        return self.field
