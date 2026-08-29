
import os
import numpy as np
import zipfile
import fidimag.common.helper as helper
from fidimag.common.integrators import CvodeSolver, CvodeSolver_OpenMP, \
    ErkSolver, StepIntegrator, ScipyIntegrator


class DriverBase:
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
        Set the Gilbert damping of the system, as a uniform or spatially
        dependent scalar field.

        Parameters
        ----------
        value
            The damping, specified in any of these ways:

            - A float from 0 to 1, for a uniform damping across the whole
              sample.

            - A function of the spatial coordinates returning values from 0
              to 1. For example, a damping that increases linearly in the x
              direction::

                  def alpha_profile(r):
                      if r[0] <= 10:
                          return r[0] / 10.
                      else:
                          return 0

            - An array with n values from 0 to 1, in the same order as the
              mesh coordinates array.

            - An array previously saved to a numpy file, loaded with
              ``numpy.load``.
        """
        self._alpha[:] = helper.init_scalar(value, self.mesh)

    alpha = property(get_alpha, set_alpha)

    def set_integrator(self, integrator, use_jac):
        """
        Choose the method used to evolve the equation in time.

        The right hand side comes from the `sundials_rhs` of the driver in
        use (LLG, LLG_STT, and so on); this only selects what integrates it.

        Parameters
        ----------
        integrator
            One of:

            - ``'sundials'`` (default): CVODE with backward differentiation
              formulae and a Newton iteration, solved with restarted GMRES.
              The implicit choice, and the one to keep on a fine mesh: the
              stiffness of the LLG equation comes from the exchange
              interaction, whose fastest timescale grows as ``1 / dx ** 2``.
            - ``'sundials_diag'``: the same, with CVODE's diagonal
              approximate Jacobian instead of GMRES.
            - ``'sundials_adams'``: CVODE with Adams-Moulton and a fixed
              point iteration. The non-stiff arm of the same solver, with no
              linear solve at all.
            - ``'dopri5'``: explicit Runge-Kutta 5(4) through ARKODE, using
              the Dormand and Prince tableau. This is the same method OOMMF
              calls ``rkf54m``, so the two codes can be compared directly.
              Note that OOMMF's own default, ``rkf54``, is RK5(4)7FC, another
              member of that family rather than the Fehlberg tableau its name
              suggests.
            - ``'rkf45'``: explicit Runge-Kutta through ARKODE with the
              genuine Fehlberg 5(4) tableau.
            - ``'euler'``, ``'rk4'``: fixed step, mostly for debugging.
            - ``'scipy'``: SciPy's integrator.
            - the ``_openmp`` variants of the sundials ones, which use the
              OpenMP N_Vector.
        use_jac
            Supply the analytical Jacobian-times-vector product to CVODE.
            Only used by the ``'sundials'`` and ``'sundials_openmp'`` options.

        Notes
        -----
        On muMAG standard problem 4, all four of the SUNDIALS options agree
        to within 2e-8 of each other and put the crossing of ``<m_x> = 0`` at
        the same place, at rtol = atol = 1e-10. Their cost differs: taking
        the BDF default as 1, Adams ran in 0.66 of the time and the two
        explicit methods in 0.47. The explicit ones evaluate the right hand
        side more often, but a step costs only its stages, with no Newton
        iteration and no Krylov solve.

        That ordering is not general. It reflects a problem whose dynamics
        are worth resolving anyway; on a finer mesh, or when relaxing to
        equilibrium rather than following the dynamics, the implicit methods
        recover their advantage. For equilibrium, prefer the minimisers.
        """

        if integrator == "sundials" and use_jac:
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs,
                                          self.sundials_jtimes)
        elif integrator == "sundials_diag":
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs,
                                          linear_solver="diag")
        elif integrator == "sundials":
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs)
        elif integrator == "sundials_adams":
            # Adams-Moulton with a fixed point iteration: the non-stiff arm of
            # CVODE, with no linear solve at all
            self.integrator = CvodeSolver(self.spin, self.sundials_rhs,
                                          lmm="adams")
        elif integrator in ("dopri5", "rkf45", "erk"):
            # Explicit Runge-Kutta, through ARKODE. `dopri5` is the tableau
            # OOMMF calls rkf54m, so the two codes can be compared directly;
            # `rkf45` is the genuine Fehlberg one
            table = {"dopri5": "ARKODE_DORMAND_PRINCE_7_4_5",
                     "erk": "ARKODE_DORMAND_PRINCE_7_4_5",
                     "rkf45": "ARKODE_FEHLBERG_6_4_5"}[integrator]
            self.integrator = ErkSolver(self.spin, self.sundials_rhs,
                                        table=table)
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
            raise NotImplementedError(
                "integrator is {}, should be one of sundials, "
                "sundials_diag, sundials_adams, dopri5, rkf45, euler, "
                "rk4, scipy, or their _openmp variants".format(integrator))

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
