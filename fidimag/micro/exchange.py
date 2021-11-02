import fidimag.extensions.micro_clib as micro_clib
from .energy import Energy
#from constant import mu_0


class UniformExchange(Energy):

    """
    UniformExchange(A, name='UniformExchange')

    Compute the exchange field in micromagnetics.
    
    Inputs:
        A: float
            A is the exchange stiffness constant measured in 
            Joules / Meter (J / M)
    
    """

    def __init__(self, A, name='UniformExchange'):
        self.A = A
        self.name = name
        self.jac = True

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        micro_clib.compute_exchange_field_micro(m,
                                                self.field,
                                                self.energy,
                                                self.Ms_inv,
                                                self.A,
                                                self.dx,
                                                self.dy,
                                                self.dz,
                                                self.n,
                                                self.neighbours
                                                )

        return self.field


class UniformAnisotropicExchange(Energy):
    r"""
    UniformAnisotropicExchange(gam, name='UniformAnisotropicExchange')

    Compute the anisotropic exchange field in micromagnetics, defined as

                              ---    d^2 (m_i)
        H_ex =  2 * gam       \      ----
                -------  *    /__    d i^2
                mu0 * Ms    i=x,y,z
    Inputs:
        gam: float
             gam is the anisotropic exchange constant measured in 
             Joules / Meter (J / m)
    
    """

    def __init__(self, gam, name='UniformAnisotropicExchange'):
        self.gam = gam
        self.name = name
        self.jac = True

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        micro_clib.compute_exchange_field_anisotropic_micro(m,
                                                            self.field,
                                                            self.energy,
                                                            self.Ms_inv,
                                                            self.gam,
                                                            self.dx,
                                                            self.dy,
                                                            self.dz,
                                                            self.n,
                                                            self.neighbours
                                                            )

        return self.field
