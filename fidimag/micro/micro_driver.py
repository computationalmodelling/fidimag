import numpy as np
from fidimag.common.integrators import CvodeSolver, CvodeSolver_OpenMP, StepIntegrator
from fidimag.common.fileio import DataSaver, DataReader
from fidimag.common.save_vtk import SaveVTK
import time
import os


class MicroDriver(object):
    """

    A class with shared methods and properties for different drivers to solve
    the Landau-Lifshitz-Gilbert equation

    Variables that are proper of the driver class:

        * Ms_const
        * t
        * spin_last
        * dm_dt
        * integrator_tolerances_set
        * step
        * n (from mesh)
        * n_nonzero (from mesh and set_Ms method in Simulation)
        * gamma
        * do_precession
        * default_c (correction factor to keep the magnetisation normalised
                     during the LLG equation integration. See the
                     fidimag/atomistic/lib/llg.c file for more details)

                     TODO: Check default_c units in the micromagnetic
                           context

    """

    def __init__(self, mesh, spin, Ms, field, alpha, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False
                 ):

        # These are (ideally) references to arrays taken from the Simulation
        # class. Variables with underscore are arrays changed by a property in
        # the simulation class
        self.mesh = mesh
        self.spin = spin
        self._Ms = Ms

        # Only for LLG STT: (??)
        self.Ms_const = 0

        self.field = field
        self._alpha = alpha
        self._pins = pins
        self.interactions = interactions
        # Strings are not referenced, this is a copy:
        self.name = name

        # The following are proper of the driver class: -----------------------
        # See also the set_default_options() function

        self.t = 0
        self.spin_last = np.ones(3 * self.mesh.n, dtype=np.float)
        self.dm_dt = np.zeros(3 * self.mesh.n, dtype=np.float)
        self.integrator_tolerances_set = False
        self.step = 0

        self.n = self.mesh.n
        self.n_nonzero = self.mesh.n  # number of spins that are not zero
                                      # We check this in the set_Ms function

        # To save VTK files:
        self.vtk = SaveVTK(self.mesh, name=name)

        self.set_default_options()

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
        else:
            raise NotImplemented("integrator must be sundials, euler or rk4")

        # ---------------------------------------------------------------------

        # Initialise the table for the data file with the simulation
        # information:

        self.data_saver = data_saver

        # This should not be necessary:
        # self.data_saver.entities['skx_num'] = {
        #     'unit': '<>',
        #     'get': lambda sim: sim.skyrmion_number(),
        #     'header': 'skx_num'}

        self.data_saver.entities['rhs_evals'] = {
            'unit': '<>',
            'get': lambda sim: self.integrator.rhs_evals(),
            'header': 'rhs_evals'}

        self.data_saver.entities['real_time'] = {
            'unit': '<s>',
            'get': lambda _: time.time(),  # seconds since epoch
            'header': 'real_time'}

        self.data_saver.update_entity_order()

        # ---------------------------------------------------------------------

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
        self._Ms[:] = Ms

        self.gamma = gamma
        self.do_precession = True

    def sundials_rhs(self, t, y, ydot):
        """
        Defined in the corresponding driver class
        """
        pass

    def compute_effective_field(self, t):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the micromagnetic Energy class
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

    def stat(self):
        return self.integrator.stat()

    def reset_integrator(self, t=0):
        """
        Reset the CVODE integrator and set the simulation time to `t`
        The simulation step is reset to zero
        """
        self.integrator.reset(self.spin, t)
        self.t = t  # also reinitialise the simulation time and step
        self.step = 0

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

        self.spin_last[:] = self.spin[:]

        flag = self.integrator.run_until(t)
        if flag < 0:
            raise Exception("Run run_until failed!!!")

        self.spin[:] = self.integrator.y[:]
        self.t = t
        self.step += 1

        self.compute_effective_field(t)  # update fields before saving data
        self.data_saver.save()

    def relax(self, dt=1e-11, stopping_dmdt=0.01, max_steps=1000,
              save_m_steps=100, save_vtk_steps=100):
        """

        Evolve the system until meeting the dmdt < stopping_dmdt criteria. We
        can specify to save VTK and NPY files with the magnetisation vector
        field, every certain number of the integrator steps

        TODO: Check what dt is exactly doing

        """

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

    # -------------------------------------------------------------------------
    # Save functions ----------------------------------------------------------
    # -------------------------------------------------------------------------

    def save_vtk(self):
        """
        Save a VTK file with the magnetisation vector field (vector data)
        and the saturation magnetisation values (scalar data) as
        cell data
        """
        self.vtk.save_vtk(self.spin.reshape(-1, 3), self._Ms, step=self.step)

    def save_m(self):
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
