from __future__ import division
from __future__ import print_function

import os
import numpy as np
import zipfile
import fidimag.common.helper as helper
from fidimag.common.integrators import CvodeSolver, CvodeSolver_OpenMP, \
    StepIntegrator, ScipyIntegrator


class DriverBase(object):
    """
    Common methods for the micromagnetic and atomistic driver classes
    """

    def __init__(self):
        pass

    # -------------------------------------------------------------------------

    def initiate_variables(self, n_spins):
        """
        Common variables for both micro and atomistic drivers
        """
        self._alpha = np.zeros(n_spins, dtype=np.float64)

        self.t = 0
        self.spin_last = np.ones(3 * n_spins, dtype=np.float64)
        self.dm_dt = np.zeros(3 * n_spins, dtype=np.float64)
        self.integrator_tolerances_set = False
        self.step = 0

    def get_alpha(self):
        """
        Returns the array with the spatially dependent Gilbert damping
        per mesh/lattice site
        """
        return self._alpha

    def set_alpha(self, value):
        """

        Set the Gilbert damping of the system as a uniform or spatially
        dependent scalar field

        ARGUMENTS:

        value     :: * For a uniform damping across the whole sample, just
                       specify a float ranging from 0 to 1

                     * In addition, you can specify a function that returns
                       values ranging from 0 to 1, which depends on the spatial
                       coordinates. For example, a damping that increases
                       linearly in the x direction:

                        def alpha_profile(r):
                            for r[0] <= 10:
                                return r[0] / 10.
                            else:
                                return 0

                     * You can also manually specify an array with n values
                     ranging from 0 to 1 with the damping values, in the same
                     order than the mesh coordinates array.

                     * Alternatively, if you previously saved the damping
                     field array to a numpy file, you can load it using
                     numpy.load(my_array)
        """
        self._alpha[:] = helper.init_scalar(value, self.mesh)

    alpha = property(get_alpha, set_alpha)

    def set_integrator(self, integrator, use_jac):
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
            self.integrator = StepIntegrator(self.spin, self.step_rhs, step=integrator)
        elif integrator == "scipy":
            self.integrator = ScipyIntegrator(self.spin, self.step_rhs)

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
        else:
            raise NotImplemented("integrator must be sundials, euler or rk4")

    # ------------------------------------------------------------------------

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
        m0 = self.spin_last
        m1 = self.spin
        dm = (m1 - m0).reshape((3, -1))
        max_dm = np.max(np.sqrt(np.sum(dm ** 2, axis=0)))
        max_dmdt = max_dm / dt
        return max_dmdt

    def run_until(self, t):
        """
        Evolve the system with a micromagnetic driver (LLG, LLG_STT, etc.)
        until a specific time `t`, using the specified integrator.
        The integrator was specified with the right hand side of the
        driver equation

        """

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
                return
            else:
                raise ValueError("t must be >= sim.t")

        ode = self.integrator

        self.spin_last[:] = self.spin[:]

        flag = ode.run_until(t)

        if flag < 0:
            raise Exception("Run cython run_until failed!!!")

        self.spin[:] = ode.y[:]

        self.t = t
        self.step += 1

        # Update field before saving data
        self.compute_effective_field(t)
        self.data_saver.save()

    def relax(self, dt=10e-12, stopping_dmdt=0.01, max_steps=1000,
              save_m_steps=100, save_vtk_steps=100,
              printing=True):
        """
        Evolve the system until meeting the `dmdt` < `stopping_dmdt` criteria.

        The magnetisation dynamics will be checked a maximum of `max_steps`
        times at an interval of at least `dt` and compared to `stopping_dmdt`
        which is given in units of degrees per nanosecond.

        With `save_m_steps` and `save_vtk_steps` the magnetisation will be
        saved every given number of integrator steps to npy resp. vtk files.

        """
        while self.step < max_steps:
            _dt = max(dt, self.integrator.get_current_step())
            self.run_until(self.t + _dt)
            # Explanation: For very small integrator steps, there is no use in
            # checking for relaxation at every step and `dt` will provide the lower
            # boundary of what is acceptable (every 10 picoseconds per default).
            # On the other hand, when the integrator steps get larger than `dt`
            # we might as well let the integrator do its work uninterrupted.

            if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                self.save_vtk()
            if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                self.save_m()

            dmdt = self.compute_dmdt(_dt)
            if printing:
                print("#{:<4} t={:<8.3g} dt={:.3g} max_dmdt={:.3g}".format(
                    self.step,  # incremented in self.run_until (called above)
                    self.t,
                    _dt,
                    dmdt / self._dmdt_factor))
            if dmdt < stopping_dmdt * self._dmdt_factor:
                break

        if save_m_steps is not None:
            self.save_m()
        if save_vtk_steps is not None:
            self.save_vtk()

    # -------------------------------------------------------------------------
    # Save functions ----------------------------------------------------------
    # -------------------------------------------------------------------------

    def save_vtk(self):
        pass

    def save_m(self, ZIP=False):
        """
        Save the magnetisation/spin vector field as a numpy array in
        a NPY file. The files are saved in the `{name}_npys` folder, where
        `{name}` is the simulation name, with the file name `m_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """

        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)
        if ZIP:
            with zipfile.ZipFile('%s_m.zip'%self.name, 'a') as myzip:
                myzip.write(name)
            try:
                os.remove(name)
            except OSError:
                pass

    def save_skx(self):
        """
        Save the skyrmion number density (sk number per mesh site)
        as a numpy array in a NPY file.
        The files are saved in the `{name}_skx_npys` folder, where
        `{name}` is the simulation name, with the file name `skx_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """
        if not os.path.exists('%s_skx_npys' % self.name):
            os.makedirs('%s_skx_npys' % self.name)
        name = '%s_skx_npys/m_%g.npy' % (self.name, self.step)

        # The _skx_number array is defined in the SimBase class in Common
        np.save(name, self._skx_number)
