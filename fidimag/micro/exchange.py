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
