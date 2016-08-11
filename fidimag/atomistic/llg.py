from __future__ import division
from __future__ import print_function

import fidimag.extensions.clib as clib

from fidimag.common.llg_driver import LLG_Driver


class LLG(LLG_Driver):

    """

    This class is the driver to solve the Landau Lifshitz Gilbert equation
    which has the form:


          ds        -gamma
         ---- =    --------  ( s X H_eff  + a * s X ( s X H_eff ) )
          dt             2
                  ( 1 + a  )

    by using the Sundials library with CVODE.

    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.MicroDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

    """
    def set_default_options(self, gamma=1.0, mu_s=1.0, alpha=0.1):
        """
        Default option for the integrator
        Default gamma is for a free electron
        """
        self.default_c = -1.0
        self._alpha[:] = alpha

        # When we create the simulation, mu_s is set to the default value. This
        # is overriden when calling the set_ms method from the Siulation class
        # or when setting Ms directly (property)
        self._magnitude[:] = mu_s

        self.gamma = gamma
        self.do_precession = True
        self._dmdt_factor = 1.0

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        clib.compute_llg_rhs(ydot,
                             self.spin,
                             self.field,
                             self._alpha,
                             self._pins,
                             self.gamma,
                             self.n,
                             self.do_precession,
                             self.default_c)

        #ydot[:] = self.dm_dt[:]

        return 0

if __name__ == '__main__':
    pass
