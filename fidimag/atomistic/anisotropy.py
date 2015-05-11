import fidimag.extensions.clib as clib
import numpy as np
from energy import Energy
import fidimag.common.helper as helper


class Anisotropy(Energy):

    """
        compute the anisotropy field with the energy density E = - K (m.u)^2
    """

    def __init__(self, Ku, axis=(1, 0, 0), name='anis', direction=None):
        self.Ku = Ku
        self.name = name
        self.axis = axis
        self.jac = True

        if direction is not None:
            self.axis = direction

    def setup(self, mesh, spin, mu_s):
        super(Anisotropy, self).setup(mesh, spin, mu_s)

        self._Ku = helper.init_scalar(self.Ku, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, True)

    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_anisotropy(m,
                                self.field,
                                self.energy,
                                self._Ku,
                                self._axis,
                                self.nx,
                                self.ny,
                                self.nz)

        return self.field * self.mu_s_inv
