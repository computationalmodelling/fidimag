import numpy as np
from .driver_base import DriverBase

import fidimag.extensions.clib as clib
from fidimag.common.integrators import CvodeSolver, CvodeSolver_OpenMP, StepIntegrator


class LLG_Driver(DriverBase):

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

    def __init__(self, mesh, spin, magnitude, pins, 
                interactions, 
                field, 
                data_saver,
                integrator = "sundials", 
                use_jac=False):
        # Inherit from the driver class
        super(LLG_Driver, self).__init__(mesh, spin, magnitude, pins, interactions, field)

        self._alpha = np.zeros(self.n, dtype=np.float)

        self.data_saver = data_saver

        # Integrator options --------------------------------------------------

        # Here we set up the CVODE integrator from Sundials to evolve a
        # specific micromagnetic equation. The equations are specified in the
        # sundials_rhs function from any of the micromagnetic drivers in the
        # micromagnetic folder (LLG, LLG_STT, etc.)

        if integrator == "sundials" and use_jac:
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs,
                                          self.sundials_jtimes)
        elif integrator == "sundials_diag":
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs,
                                          linear_solver="diag")
        elif integrator == "sundials":
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs)
        elif integrator == "euler" or integrator == "rk4":
            self.integrator = CvodeSolver(self.spin, self.step_rhs,
                                          integrator)

        elif integrator == "sundials_openmp" and use_jac:
            self.integrator = CvodeSolver_OpenMP(self.spin, self.sundials_rhs,
                                                 self.sundials_jtimes)
        elif integrator == "sundials_diag_openmp":
            self.integrator = CvodeSolver_OpenMP(self.spin, self.sundials_rhs,
                                                 linear_solver="diag")
        elif integrator == "sundials_openmp":
            self.integrator = CvodeSolver_OpenMP(self.spin, self.sundials_rhs)
        elif integrator == "euler_openmp" or integrator == "rk4_openmp":
            self.integrator = CvodeSolver_OpenMP(self.spin, self.step_rhs,
                                                 integrator)
        elif integrator == "stochastic":
            self.integrator = None
            #we don't do anything here!
            pass
        else:
            raise NotImplemented("integrator must be sundials, euler or rk4")


        # ---------------------------------------------------------------------

        # OOMMF convention is to check if the spins have moved by ~1 degree in
        # a nanosecond in order to stop a simulation, so we set this scale for
        # dm/dt

        # ONE_DEGREE_PER_NANOSECOND:
        self._dmdt_factor = (2 * np.pi / 360) / 1e-9

        self.set_default_options()


    def set_default_options(self):
        pass

    def set_initial_value(self, t):
        self.t = t 
        if self.integrator is not None:
            self.integrator.set_initial_value(self.spin, t)
        else:
            pass

    # -------------------------------------------------------------------------

    def set_tols(self, rtol=1e-8, atol=1e-10, max_ord=None, reset=True):
        """
        Set the relative and absolute tolerances for the CVODE integrator
        """
        if max_ord is not None:
            self.integrator.set_options(rtol=rtol, atol=atol, max_ord=max_ord)
        else:
            # not all integrators have max_ord (only VODE does)
            # and we don't want to encode a default value here either
            self.integrator.set_options(rtol=rtol, atol=atol)
        if reset:
            self.reset_integrator(self.t)

    def stat(self):
        return self.integrator.stat()

    def set_default_options(self):
        pass

    def reset_integrator(self, t=0):
        """
        Reset the CVODE integrator and set the simulation time to `t`
        The simulation step is reset to zero
        """
        self.integrator.reset(self.spin, t)
        self.t = t  # also reinitialise the simulation time and step
        self.step = 0

    def set_tols(self, rtol=1e-8, atol=1e-10):
        """
        Set the relative and absolute tolerances for the CVODE integrator
        """
        self.integrator.set_options(rtol, atol)

    def compute_effective_field(self, t):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the corresponding Energy class
        """

        # self.spin[:] = y[:]

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t)

    def compute_effective_field_jac(self, t, spin):
        self.field[:] = 0
        for obj in self.interactions:
            if obj.jac:
                self.field += obj.compute_field(t, spin=spin)

    def compute_dmdt(self, dt):
        m0 = self.previous_spin
        m1 = self.spin
        dm = (m1 - m0).reshape((3, -1))
        max_dm = np.max(np.sqrt(np.sum(dm ** 2, axis=0)))
        max_dmdt = max_dm / dt
        return max_dmdt

    def cvode_dt(self):
        return self.integrator.get_current_step()

    def next_step(self, time=None):
        dt = 0
        if time is None:
            dt = self.cvode_dt()
            self.t += dt
        else:
            self.t = time

        self.previous_spin[:] = self.spin[:]
        flag = self.integrator.run_until(self.t)

        if flag < 0:
            raise Exception("Run cython run_until failed!!!")

        self.step += 1
        
        self.spin[:] = self.integrator.y[:]

if __name__ == '__main__':
    pass
