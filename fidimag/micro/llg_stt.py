from __future__ import division

import fidimag.extensions.common_clib as clib
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

    The current can be set in several ways.
    If the current is constant, it can be set by

    driver.jx = 1.0
    driver.jy = 0.0
    driver.jz = 0.0

    The current can also be set using a function by setting the parameter:

        driver.jx_function = func_x
        driver.jy_function = func_y
        driver.jz_function = func_z

    In order to use this, after creating the sim object, you must set the function parameters
    with functions as:

    sim.driver.jx_func = myfunc_x
    sim.driver.jy_func = myfunc_y
    sim.driver.jz_func = myfunc_z

    The function definition to set the current must be of the form:
        def myfunc_x(pos):
            x, y, z = pos
            jx = # some function of x, y, and z.
            return jx

    Or, if the current is time dependent:
        def myfunc_x(pos, t):
            x, y, z = pos
            jx = # some function of x, y, z and t.
            return jx

    Note that if not set, then *no current will be applied in that direction*.

    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.MicroDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

        """

    def __init__(self, mesh, spin, Ms, Ms_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False
                 ):

        # Inherit from the driver class
        super(LLG_STT, self).__init__(mesh, spin, Ms, Ms_inv, field,
                                      pins, interactions, name,
                                      data_saver,
                                      integrator=integrator,
                                      use_jac=use_jac
                                      )

        self.field_stt = np.zeros(3 * self.n)

        self._jx = np.zeros(self.n, dtype=np.float)
        self._jy = np.zeros(self.n, dtype=np.float)
        self._jz = np.zeros(self.n, dtype=np.float)

        self.p = 0.5
        self.beta = 0
        self.jx_function = None
        self.jy_function = None
        self.jz_function = None

        # FIXME: change the u0 to spatial
        self.u0 = const.g_e * const.mu_B / (2 * const.c_e)

    def get_jx(self):
        return self._jx

    def set_jx(self, value, *args):
        self._jx[:] = helper.init_scalar(value, self.mesh, *args)

    jx = property(get_jx, set_jx)

    def get_jy(self):
        return self._jy

    def set_jy(self, value, *args):
        self._jy[:] = helper.init_scalar(value, self.mesh, *args)

    jy = property(get_jy, set_jy)

    def get_jz(self):
        return self._jz

    def set_jz(self, value, *args):
        self._jz[:] = helper.init_scalar(value, self.mesh, *args)

    jz = property(get_jz, set_jz)

    def sundials_rhs(self, t, y, ydot):
        self.t = t
        # already synchronized when call this funciton
        # self.spin[:]=y[:]
        self.compute_effective_field(t)
        if self.jx_function:
            self.set_jx(self.jx_function, t)
        if self.jy_function:
            self.set_jy(self.jy_function, t)
        if self.jz_function:
            self.set_jz(self.jz_function, t)

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
