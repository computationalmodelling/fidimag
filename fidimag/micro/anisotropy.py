import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from fidimag.common.constant import mu_0
from .energy import Energy
import fidimag.common.helper as helper


class UniaxialAnisotropy(Energy):

    """
        compute the anisotropy field with the energy density E = K[1- (m.u)^2]
    """

    def __init__(self, Ku, axis=(1, 0, 0), name='anis'):
        self.Ku = Ku
        self.name = name
        self.jac = True
        self.axis = axis

    def setup(self, mesh, spin, Ms):
        super(UniaxialAnisotropy, self).setup(mesh, spin, Ms)

        self._Ku = helper.init_scalar(self.Ku, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, True)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        micro_clib.compute_anisotropy_micro(m,
                                            self.field,
                                            self.energy,
                                            self.Ms_inv,
                                            self._Ku,
                                            self._axis,
                                            self.nx,
                                            self.ny,
                                            self.nz)
        return self.field
