import fidimag.extensions.clib as clib
import numpy as np
from energy import Energy


class DemagFull(Energy):
    """
    Calculation of the demag field using a brute-force approach,
    summing the dipolar contributions of the whole system for
    every lattice point.
    Since we can obtain the field directly, the energy is calculated
    in this process, thus we inherit from the Energy class
    as in the Exchange field, to calculate the total energy summing up
    the energy density.
    """

    def __init__(self, name='demag_full'):
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s):
        super(DemagFull, self).setup(mesh, spin, mu_s)

        unit_length = mesh.unit_length
        self.mu_s_scale = np.zeros(mesh.n, dtype=np.float)

        # note that the 1e-7 comes from \frac{\mu_0}{4\pi}
        self.scale = 1e-7 / unit_length ** 3

        # could be wrong, needs carefully tests!!!
        self.mu_s_scale = self.mu_s * self.scale

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_demag_full(m, self.field, self.mesh.coordinates,
                                self.energy, self.mu_s, self.n)

        return self.field * self.mu_s_scale

    # def compute_energy(self):
    #     energy = self.demag.compute_energy(
    #         self.spin, self.mu_s_scale, self.field)
    #     return energy / self.scale
