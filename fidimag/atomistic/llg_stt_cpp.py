from __future__ import division

import types

import fidimag.extensions.clib as clib
import numpy as np

from fidimag.common.llg_driver import LLG_Driver
import fidimag.common.helper as helper


class LLG_STT_CPP(LLG_Driver):

    """

    This class is the driver to solve the Landau Lifshitz Gilbert equation
    with a Current Perpendicular to the Plane, which follows the formalism
    of Spin Transfer Torque. The equation is given by:


          ds        -gamma
         ---- =    --------  ( s X H_eff  + a * s X ( s X H_eff ) ) + .... ADD
          dt             2
                  ( 1 + a  )


    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.MicroDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

    """

    def __init__(self, mesh, spin, magnitude, pins, 
                interactions, 
                field, 
                data_saver, integrator = "sundials", 
                use_jac=False):
    
        # Inherit from the driver class
        super(LLG_STT_CPP, self).__init__(mesh, spin, magnitude, pins, 
                                        interactions, field, 
                                        data_saver,
                                        integrator = integrator, 
                                        use_jac=False)

        self._p = np.zeros(3 * self.n, dtype=np.float)
        self._a_J = np.zeros(self.n, dtype=np.float)

        self.a_J = 1
        self.beta = 0
        self.J_time_fun = None

    def get_p(self):
        return self._p

    def set_p(self, value):
        self._p[:] = helper.init_vector(value, self.mesh)
    p = property(get_p, set_p)

    def get_a_J(self):
        return self._a_J

    def set_a_J(self, value):
        self._a_J[:] = helper.init_scalar(value, self.mesh)

    a_J = property(get_a_J, set_a_J)

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        if self.J_time_fun is not None:
            clib.compute_llg_stt_cpp(ydot,
                                     self.spin,
                                     self.field,
                                     self._p,
                                     self._alpha,
                                     self._pins,
                                     self._a_J*self.J_time_fun(t),
                                     self.beta,
                                     self.gamma,
                                     self.n)
        else:
            clib.compute_llg_stt_cpp(ydot,
                                     self.spin,
                                     self.field,
                                     self._p,
                                     self._alpha,
                                     self._pins,
                                     self._a_J,
                                     self.beta,
                                     self.gamma,
                                     self.n)
