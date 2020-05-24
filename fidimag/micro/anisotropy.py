import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from fidimag.common.constant import mu_0
from .energy import Energy
import fidimag.common.helper as helper


class UniaxialAnisotropy(Energy):

    """
        compute the anisotropy field with the energy density E = K[1- (m.u)^2]
    """

    def __init__(self, Ku, axis=(1, 0, 0), name='Anisotropy'):
        self.Ku = Ku
        self.name = name
        self.jac = True
        self.axis = axis

    def setup(self, mesh, spin, Ms, Ms_inv):
        super(UniaxialAnisotropy, self).setup(mesh, spin, Ms, Ms_inv)

        self._Ku = helper.init_scalar(self.Ku, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, 3, norm=True)

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


class UniaxialAnisotropy4(Energy):

    """
    4th Order Uniaxial Anisotropy
    """

    def __init__(self, K1, K2, axis=(1, 0, 0), name='UAnisotropy4'):
        self.K1 = K1
        self.K2 = K2
        self.name = name
        self.jac = True
        self.axis = axis

    def setup(self, mesh, spin, Ms, Ms_inv):
        super(UniaxialAnisotropy4, self).setup(mesh, spin, Ms, Ms_inv)

        self._K1 = helper.init_scalar(self.K1, self.mesh)
        self._K2 = helper.init_scalar(self.K2, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, 3, norm=True)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        micro_clib.compute_anisotropy4_micro(m,
                                             self.field,
                                             self.energy,
                                             self.Ms_inv,
                                             self._K1,
                                             self._K2,
                                             self._axis,
                                             self.nx,
                                             self.ny,
                                             self.nz)





        return self.field

