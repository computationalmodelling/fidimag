import fidimag.extensions.micro_clib as micro_clib
from .energy import Energy
#from constant import mu_0


class ExchangeRKKY(Energy):

    """
    RKKYExchange(sigma, Delta, z_bottom=0, z_top = -1, name='RKKYExchange')

    Compute the RKKY-style exchange interaction defined by

    E = sigma/Delta *(1-m_t * m_b)

    where E is the energy density, sigma is the surface exchange coefficient between the two surfaces, 
    Delata is the space thickness. m_t and m_b are the unit vectors of the top and bottom layers, respectively.
    
    Inputs:
        sigma: float
            sigma is the surface exchange stiffness constant.
        Delta: float
            Delta is the the space thickness in Meter.
        z_bottom: int
            z_bottom is the index of the bottom layer
        z_top: int
            z_top is the index of the top layer.
    
    """

    def __init__(self, sigma, Delta=1e-9, z_bottom=0, z_top = -1, name='RKKYExchange'):
        self.sigma = sigma/Delta
        self.name = name
        self.z_bottom = z_bottom
        self.z_top = z_top
        self.jac = True

    def setup(self, mesh, spin, Ms, Ms_inv):
        super(ExchangeRKKY, self).setup(mesh, spin, Ms, Ms_inv)
        if self.z_top<0:
            self.z_top+= self.nz
 

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        micro_clib.compute_exchange_field_micro_rkky(m,
                                                self.field,
                                                self.energy,
                                                self.Ms_inv,
                                                self.sigma,
                                                self.nx,
                                                self.ny,
                                                self.nz,
                                                self.z_bottom,
                                                self.z_top
                                                )

        return self.field
