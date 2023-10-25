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

                __     ^     ->  4
        w = -Kc \   (  a_i . m  )
                /_
                i

    with `a_i` the three perpendicular anisotropy axes. The three axes
    are specified using the `cubicAxes` argument. Notice that in this
    definition, `Kc > 0` means `a_i` are easy axes.
    """

    def __init__(self, Kc, cubicAxes=(1., 0., 0., 0., 1., 0., 0., 0., 1.),
                 name='CubicAnisotropy'):
        self.Kc = Kc
        self.name = name
        self.jac = True
        self.cubicAxes = cubicAxes

    @property
    def cubicAxes(self):
        return self._cubicAxes

    @cubicAxes.setter
    def cubicAxes(self, axes):
        self._cubicAxes = np.array(axes).astype(np.float64)
        if self._cubicAxes.size != 9:
            raise ValueError('Sequence with axes coordinates requires 9 input values')

        self._cubicAxes.shape = (3, 3)

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(CubicAnisotropy, self).setup(mesh, spin, mu_s, mu_s_inv)
        self._Kc = helper.init_scalar(self.Kc, self.mesh)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin.reshape(-1, 3)
        else:
            m = self.spin.reshape(-1, 3)

        f_vec = self.field.reshape(-1, 3)
        e_rs = self.energy
        f_vec[:] = 0.0
        e_rs[:] = 0.0

        factor = (4. * self._Kc * self.mu_s_inv)[:, None]  # (n, 1) array
        for ax in self._cubicAxes:
            # Every row is m . ax . Make a column vector: (n, 1)
            m_dot_ax = np.array(np.einsum('ij,j->i', m, ax), ndmin=2).reshape(-1, 1)
            # factor * m_dot_ax generate a column vector. Every row is multiplied by ax
            f_vec[:] += (factor * (m_dot_ax ** 3)) * ax
            e_rs[:] += np.einsum('i,ij->i', -self._Kc, m_dot_ax ** 4)

        # clib.compute_anisotropy_cubic(m,
        #                               self.field,
        #                               self.mu_s_inv,
        #                               self.energy,
        #                               self._Kc,
        #                               self.n)

        return self.field
