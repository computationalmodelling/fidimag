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
    r"""
    Compute the cubic anisotropy field using 2 specified perpendicular axes
    directions (the 3rd is set perpendicular). The energy density reads::

                __         ^     ->  4
        w = -Kc \   m_i (  a_i . m  )
                /_
                i

    with a_i the three perpendicular anisotropy axes. Only two of the axes
    are specified using the `axisA` and `axisB` arguments.
    """

    def __init__(self, Kc, axisA=(1., 0, 0), axisB=(0, 1., 0),
                 name='CubicAnisotropy'):
        self.Kc = Kc
        self.name = name
        self.jac = True

        self.axisA = np.array(axisA).astype(np.float64)
        self.axisA /= np.linalg.norm(self.axisA)
        self.axisB = np.array(axisB).astype(np.float64)
        self.axisB /= np.linalg.norm(self.axisB)

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(CubicAnisotropy, self).setup(mesh, spin, mu_s, mu_s_inv)
        self._Kc = helper.init_scalar(self.Kc, self.mesh)
        self.axisC = np.cross(self.axisA, self.axisB)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin.reshape(-1, 3)
        else:
            m = self.spin.reshape(-1, 3)

        f_vec = self.field.reshape(-1, 3)
        factor = (4. * self._Kc * self.mu_s_inv)[:, np.newaxis]
        f_vec[:] = factor * (np.einsum('ij,j->i', m, self.axisA) ** 3)[:, np.newaxis] * self.axisA
        f_vec[:] += factor * (np.einsum('ij,j->i', m, self.axisB) ** 3)[:, np.newaxis] * self.axisB
        f_vec[:] += factor * (np.einsum('ij,j->i', m, self.axisC) ** 3)[:, np.newaxis] * self.axisC

        e_rs = self.energy.reshape(-1, 3)
        e_rs[:] = -self._Kc * (np.einsum('ij,j->i', m, self.axisA) ** 4)[:, np.newaxis]
        e_rs[:] += -self._Kc * (np.einsum('ij,j->i', m, self.axisB) ** 4)[:, np.newaxis]
        e_rs[:] += -self._Kc * (np.einsum('ij,j->i', m, self.axisC) ** 4)[:, np.newaxis]

        # clib.compute_anisotropy_cubic(m,
        #                               self.field,
        #                               self.mu_s_inv,
        #                               self.energy,
        #                               self._Kc,
        #                               self.n)

        return self.field
