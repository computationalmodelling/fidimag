from __future__ import division

import fidimag.extensions.clib as clib
import numpy as np

from .micro_driver import MicroDriver
import fidimag.common.helper as helper
import fidimag.common.constant as const


class LLG_STT(MicroDriver):

    """

    This class is the driver to solve the Landau Lifshitz Gilbert equation
    with a Spin Transfer Torque term, which has the form:


      dm        -gamma
     ---- =    --------  ( m X H_eff  + a * m X ( m x H_eff ) + ... )
      dt             2
              ( 1 + a  )

    by using the Sundials library with CVODE.

    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.MicroDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

        """

    def __init__(self, mesh, spin, Ms, field, alpha, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False
                 ):

        # Inherit from the driver class
        super(LLG_STT, self).__init__(mesh, spin, Ms, field,
                                      alpha, pins, interactions, name,
                                      data_saver,
                                      integrator='sundials',
                                      use_jac=False
                                      )

        self.field_stt = np.zeros(3 * self.n)

        self._jx = np.zeros(self.n, dtype=np.float)
        self._jy = np.zeros(self.n, dtype=np.float)
        self._jz = np.zeros(self.n, dtype=np.float)

        self.p = 0.5
        self.beta = 0
        self.update_j_fun = None

        # FIXME: change the u0 to spatial
        self.u0 = const.g_e * const.mu_B / (2 * const.c_e)

    def get_jx(self):
        return self._jx

    def set_jx(self, value):
        self._jx[:] = helper.init_scalar(value, self.mesh)

    jx = property(get_jx, set_jx)

    def get_jy(self):
        return self._jy

    def set_jy(self, value):
        self._jy[:] = helper.init_scalar(value, self.mesh)

    jy = property(get_jy, set_jy)

    def get_jz(self):
        return self._jz

    def set_jz(self, value):
        self._jz[:] = helper.init_scalar(value, self.mesh)

    jz = property(get_jz, set_jz)

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        if self.update_j_fun is not None:
            clib.compute_stt_field(self.spin,
                                   self.field_stt,
                                   self._jx * self.update_j_fun(t),
                                   self._jy * self.update_j_fun(t),
                                   self._jz * self.update_j_fun(t),
                                   self.mesh.dx * self.mesh.unit_length,
                                   self.mesh.dy * self.mesh.unit_length,
                                   self.mesh.dz * self.mesh.unit_length,
                                   self.mesh.neighbours,
                                   self.n
                                   )
        else:
            clib.compute_stt_field(self.spin,
                                   self.field_stt,
                                   self._jx,
                                   self._jy,
                                   self._jz,
                                   self.mesh.dx * self.mesh.unit_length,
                                   self.mesh.dy * self.mesh.unit_length,
                                   self.mesh.dz * self.mesh.unit_length,
                                   self.mesh.neighbours,
                                   self.n
                                   )

        clib.compute_llg_stt_rhs(ydot,
                                 self.spin,
                                 self.field,
                                 self.field_stt,
                                 self._alpha,
                                 self.beta,
                                 self.u0 * self.p / self.Ms_const,
                                 self.gamma,
                                 self.n)
