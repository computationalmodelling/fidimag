from __future__ import division
from __future__ import print_function
import fidimag.extensions.clib as clib

from micro_driver import MicroDriver


class LLG(MicroDriver):
    """

    This class is the driver to solve the Landau Lifshitz Gilbert equation
    which has the form:


          dm        -gamma
         ---- =    --------  ( m X H_eff  + a * m X ( m x H_eff ) )
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
        super(LLG, self).__init__(mesh, spin, Ms, field,
                                  alpha, pins, interactions, name,
                                  data_saver,
                                  integrator='sundials',
                                  use_jac=False
                                  )

    def sundials_rhs(self, t, y, ydot):
        """
        The right hand side of the LLG equation
        """

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
                             self.mesh.n,
                             self.do_precession,
                             self.default_c)

        # ydot[:] = self.dm_dt[:]

        return 0

    def sundials_jtimes(self, mp, Jmp, t, m, fy):
        # From the micro_driver class:
        self.compute_effective_field_jac(t, mp)

        clib.compute_llg_jtimes(Jmp,
                                m, fy,
                                mp, self.field,
                                self.alpha,
                                self._pins,
                                self.gamma,
                                self.n,
                                self.do_precession,
                                self.default_c)
        return 0

    def step_rhs(self, t, y):
        self.spin[:] = y[:]
        self.t = t

        # From the micro_driver class:
        self.compute_effective_field(t)

        clib.compute_llg_rhs(self.dm_dt,
                             self.spin,
                             self.field,
                             self.alpha,
                             self._pins,
                             self.gamma,
                             self.n,
                             self.do_precession,
                             self.default_c)
        return self.dm_dt

if __name__ == '__main__':
    pass
