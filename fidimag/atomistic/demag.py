import fidimag.extensions.dipolar as clib
import numpy as np
from .energy import Energy
import fidimag
from fidimag.atomistic.energy import Energy
import fidimag.extensions.fmm as fmm
import time
import sys


class Demag(Energy):
    r"""

    Energy class for the demagnetising field (a.k.a. dipolar interactions,
    stray field), *only for Cuboid meshes* (i.e. a square lattice in the
    discrete spin model), since this class uses the OOMMF's FFT code to
    simplify the field calculations.

    The field has the expression::

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

    def __init__(self, calc_every=0, name='Demag'):
        self.calc_every = calc_every
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super().setup(mesh, spin, mu_s, mu_s_inv)

        # Ryan Pepper 04/04/2019
        # We *do not* need to scale by mesh.unit_length**3 here!
        # This is because in the base energy class, dx, dy and dz
        # are scaled.
        self.scale = 1e-7
        self.mu_s_scale = mu_s * self.scale

        self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                   self.nx, self.ny, self.nz,
                                   tensor_type='dipolar')
        if not self.calc_every:
            self.compute_field = self.compute_field_every
        else:
            self.count = 0
            self.compute_field = self.compute_field_periodically

    def compute_field_every(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        self.demag.compute_field(m, self.mu_s_scale, self.field)
        return self.field

    def compute_field_periodically(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        if not (self.count % self.calc_every == 0):
            self.count += 1
            return self.field
        else:
            print(self.count)
            self.count += 1
            self.demag.compute_field(m, self.mu_s_scale, self.field)
            return self.field


    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin, self.mu_s_scale, field)
        return field

    def compute_energy(self):
        # The C routine takes the field as an input, so it has to be the
        # field at the current spins. The drivers also read `field` straight
        # after asking for the energy
        self.compute_field()
        self.demag.compute_energy(
            self.spin, self.mu_s_scale, self.field, self.energy)

        # The C routine fills `energy` per site and returns their sum; the
        # scaling is undone here and the total kept, as the base class does
        self.energy /= self.scale
        self.total_energy = np.sum(self.energy)

        return self.total_energy


class DemagFMM(Energy):
    def __init__(self, order, ncrit, theta, name="DemagFMM", type='fmm'):
        self.type = type
        if self.type == 'fmm':
            self._type = 0
        elif self.type == 'bh':
            self._type = 1

        self.name = name
        assert order > 0, "Order must be 1 or higher"
        assert order < 11, "Order bust be < 11"
        self.order = order
        assert ncrit >= 2, "ncrit must be greater than 1."
        self.ncrit = ncrit
        assert theta >= 0.0, "theta must be >= 0.0"
        self.theta = theta

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super().setup(mesh, spin, mu_s, mu_s_inv)
        self.coords = mesh.coordinates * mesh.unit_length
        self.m_temp = np.zeros_like(spin)
        self.fmm = fmm.FMM(self.n, self.ncrit, self.theta,
                           self.order,
                           self.coords,
                           self.m_temp,
                           self.mu_s, self._type
                           )

    def compute_field(self, t=0, spin=None):
        self.m_temp[:] = spin if spin is not None else self.spin
        self.fmm.compute_field(self.field)
        self.field *= -1e-7
        return self.field

    def compute_energy(self):
        # Energy is computed lazily (only when requested), mirroring the FFT
        # Demag class, so that compute_field stays free of this cost on the
        # integrator hot path.
        self.compute_field()

        m = self.spin
        # Per-site energy density, same convention as DemagFull (demag_full.c):
        # E_i = -(mu_s_i / 2) * m_i . H_i , the 1/2 avoiding double counting.
        # self.field is already the scaled demag field.
        self.energy[:] = -0.5 * self.mu_s * (self.field[0::3] * m[0::3]
                                             + self.field[1::3] * m[1::3]
                                             + self.field[2::3] * m[2::3])

        self.total_energy = np.sum(self.energy)
        return self.total_energy
