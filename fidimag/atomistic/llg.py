from __future__ import division
from __future__ import print_function

import fidimag.extensions.clib as clib

from .atomistic_driver import AtomisticDriver


class LLG(AtomisticDriver):

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

    def __init__(self, mesh, spin, mu_s, mu_s_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac
                 ):

        # Inherit from the driver class
        super(LLG, self).__init__(mesh, spin, mu_s, mu_s_inv, field,
                                  pins, interactions, name,
                                  data_saver,
                                  use_jac
                                  )

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
