from __future__ import division
from __future__ import print_function
import os
import time
import fidimag.extensions.clib as clib
import numpy as np

import re

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

    For now, we will copy the spin array passed to this function and after
    relaxation or evolving the LLG eq, we return the array so it is updated in
    the main Simulation class (???)

    """

    def __init__(self, spin, alpha, field,
                 integrator='sundials',
                 use_jac=False
                 ):
        """

        """

        # Inherit from the driver class
        super(LLG, self).__init__(integrator, use_jac)

    def run_until(self, t):
        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.saver.save()
            return

        self.spin_last[:] = self.spin[:]

        flag = self.integrator.run_until(t)
        if flag < 0:
            raise Exception("Run run_until failed!!!")

        self.spin[:] = self.integrator.y[:]
        self.t = t
        self.step += 1

        self.compute_effective_field(t)  # update fields before saving data
        self.saver.save()

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        clib.compute_llg_rhs(ydot,
                             self.spin,
                             self.field,
                             self.alpha,
                             self._pins,
                             self.gamma,
                             self.n,
                             self.do_precession,
                             self.default_c)

        # ydot[:] = self.dm_dt[:]

        return 0

    def sundials_jtimes(self, mp, Jmp, t, m, fy):
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

    def relax(self, dt=1e-11, stopping_dmdt=0.01, max_steps=1000,
              save_m_steps=100, save_vtk_steps=100):

        # OOMMF convention is to check if the spins have moved by
        # ~1 degree in a nanosecond in order to stop a simulation,
        # so we set this scale for dm/dt
        ONE_DEGREE_PER_NS = (2 * np.pi / 360) / 1e-9

        for i in range(0, max_steps + 1):

            cvode_dt = self.integrator.get_current_step()

            increment_dt = dt

            if cvode_dt > dt:
                increment_dt = cvode_dt

            self.run_until(self.t + increment_dt)

            if save_vtk_steps is not None:
                if i % save_vtk_steps == 0:
                    self.save_vtk()
            if save_m_steps is not None:
                if i % save_m_steps == 0:
                    self.save_m()

            dmdt = self.compute_dmdt(increment_dt)

            print(('step=%d ' +
                   'time=%0.3g ' +
                   'max_dmdt=%0.3g ' +
                   'ode_step=%0.3g') % (self.step,
                                        self.t,
                                        dmdt / ONE_DEGREE_PER_NS,
                                        cvode_dt)
                  )

            if dmdt < stopping_dmdt * ONE_DEGREE_PER_NS:
                break

        if save_m_steps is not None:
            self.save_m()

        if save_vtk_steps is not None:
            self.save_vtk()

if __name__ == '__main__':
    pass
