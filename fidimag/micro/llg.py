from __future__ import division
from __future__ import print_function

from fidimag.common.llg_driver import LLG_Driver

# We use the atomistic/lib/llg.c file to calculate the LLG equation for
# the micromagnetic case -> Move this library to Common in the future
import fidimag.extensions.clib as clib

#from .micro_driver import Driver


class LLG(LLG_Driver):
    """

    This class is the driver to solve the Landau Lifshitz Gilbert equation
    which has the form:


          dm        -gamma
         ---- =    --------  ( m X H_eff  + a * m X ( m X H_eff ) )
          dt             2
                  ( 1 + a  )

    by using the Sundials library with CVODE.

    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.MicroDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

    """

    def set_default_options(self, gamma=2.21e5, Ms=8.0e5, alpha=0.1):
        """
        Default option for the integrator
        Default gamma is for a free electron
        """
        self.default_c = 1e11
        self._alpha[:] = alpha

        # When we create the simulation, Ms is set to the default value. This
        # is overriden when calling the set_Ms method from the Siulation class
        # or when setting Ms directly (property)
        self._magnitude[:] = Ms

        self.gamma = gamma
        self.do_precession = True

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
                             self.default_c
                             )

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
                                self.default_c
                                )
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
                             self.default_c
                             )
        return self.dm_dt

if __name__ == '__main__':
    pass
