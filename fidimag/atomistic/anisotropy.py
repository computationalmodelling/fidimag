import fidimag.extensions.clib as clib
import numpy as np
from .energy import Energy
import fidimag.common.helper as helper


class Anisotropy(Energy):

    """

    This class provides an anisotropy term to the system energy, which
    is defined as

                  __         ->    ^     2
         E =  -  \    K_i  ( S_i * u_i )
                 /__
                  i

    with K_i the anisotropy constant at the i-th site and u_i the unitary
    direction of the anisotropy vector at the i-th site, thus the magnitude and
    axis can be space dependent.

    OPTIONAL ARGUMENTS --------------------------------------------------------

        Ku          :: The anisotropy constant. Can be a constant or a space
                       dependent scalar field, given as a function or array.

        axis        :: The unitary axis vector. It can be a 3 tuple, list
                       or array (uniaxial anisotropy), or a space dependent
                       vector field, gicen as a function or array.

        name        :: Interaction name

    USAGE ---------------------------------------------------------------------

    Considering a simulation object *Sim*, an uniaxial anisotropy along
    the z direction can be defined as

            K = 1 * meV
            Sim.add(Anisotropy(K, axis=(0, 0, 1)))

    For a space dependent anisotropy along the x direction, we can define a
    scalar field for the magnitude and a vector field for the anisotropy axes.

    """

    def __init__(self, Ku, axis=(1, 0, 0), name='Anisotropy'):
        self.Ku = Ku
        self.name = name
        self.axis = axis
        self.jac = True

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(Anisotropy, self).setup(mesh, spin, mu_s, mu_s_inv)

        self._Ku = helper.init_scalar(self.Ku, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, norm=True)

    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_anisotropy(m,
                                self.field,
                                self.mu_s_inv,
                                self.energy,
                                self._Ku,
                                self._axis,
                                self.n
                                )

        return self.field


class CubicAnisotropy(Energy):
    """
    Compute the Cubic Anisotropy, see documentation for detailed equations.
    """

    def __init__(self, Kc, name='CubicAnisotropy'):
        self.Kc = Kc
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(CubicAnisotropy, self).setup(mesh, spin, mu_s, mu_s_inv)
        self._Kc = helper.init_scalar(self.Kc, self.mesh)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_anisotropy_cubic(m,
                                      self.field,
                                      self.mu_s_inv,
                                      self.energy,
                                      self._Kc,
                                      self.n)

        return self.field
