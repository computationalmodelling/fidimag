import fidimag.extensions.dipolar as clib
import numpy as np
from .energy import Energy


class Demag(Energy):
    """

    Energy class for the demagnetising field (a.k.a. dipolar interactions,
    stray field), *only for Cuboid meshes* (i.e. a square lattice in the
    discrete spin model), since this class uses the OOMMF's FFT code to
    simplify the field calculations.

    The field has the expression:

                                   ^      ^         ^        ^
       ->      mu0 mu_s    __    3 r_ij ( m_j \cdot r_ij ) - m_j
       H_i =   --------   \	   -------------------------------
                 4 pi     /__              r_ij ^ 3

                        i != j

     where the numerator has unit vectors (^) and
                                                     ->    ->
     r_ij is a vector from r_i to r_j , i.e.  r_ij = r_j - r_i

     Accordingly, the energy is computed as:

                 mu_s    __   ^         ->
       E_i =  -   --    \     m_i \cdot H_i
                   2    /__

                      i=x,y,z

    """

    def __init__(self, name='Demag'):
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(Demag, self).setup(mesh, spin, mu_s, mu_s_inv)
        self.scale = 1e-7 / mesh.unit_length**3

        # could be wrong, needs carefully tests!!!
        # David Tue 19 Jun 2018: This variable is updated in the SIM class in
        # case mu_s changes
        self.mu_s_scale = mu_s * self.scale

        self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                   self.nx, self.ny, self.nz,
                                   tensor_type='dipolar')

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        self.demag.compute_field(m, self.mu_s_scale, self.field)
        return self.field

    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin, self.mu_s_scale, field)
        return field

    def compute_energy(self):

        energy = self.demag.compute_energy(
            self.spin, self.mu_s_scale, self.field, self.energy)

        return energy / self.scale
