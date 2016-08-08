from __future__ import division
from __future__ import print_function

import fidimag.extensions.clib as clib
import fidimag.extensions.cvode as cvode
from fidimag.common.save_vtk import SaveVTK
import fidimag.common.constant as const

import numpy as np
import os


class AtomisticDriver(object):
    """

    A class with shared methods and properties for different drivers to solve
    the Landau-Lifshitz-Gilbert equation

    Variables that are proper of the driver class:

        * mu_s_const
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

    """

    def __init__(self, mesh, spin, mu_s, mu_s_inv, field, alpha, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False
                 ):

        # These are (ideally) references to arrays taken from the Simulation
        # class. Variables with underscore are arrays changed by a property in
        # the simulation class
        self.mesh = mesh
        self.spin = spin
        self._mu_s = mu_s
        self._mu_s_inv = mu_s_inv

        # Only for LLG STT: (??)
        self.mu_s_const = 0

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

        if use_jac is not True:
            self.integrator = cvode.CvodeSolver(self.spin, self.sundials_rhs)
        else:
            self.integrator = cvode.CvodeSolver(
                self.spin, self.sundials_rhs, self.sundials_jtn)

        self.set_default_options()

        self.set_tols()

        # When initialising the integrator in the self.integrator call, the CVOde
        # class calls the set_initial_value function (with flag_m=0), which
        # initialises a new integrator and allocates memory in this process.
        # Now, when we set the magnetisation, we will use the same memory
        # setting this flag_m to 1, so instead of calling CVodeInit we call
        # CVodeReInit. If don't, memory is allocated in every call of set_m
        self.flag_m = 1

        # ---------------------------------------------------------------------

        self.data_saver = data_saver

        # Initialise the table for the data file with the simulation
        # information:

        # self.saver.entities['skx_num'] = {
        #     'unit': '<>',
        #     'get': lambda sim: sim.skyrmion_number(),
        #     'header': 'skx_num'}

        # self.saver.update_entity_order()

        # ---------------------------------------------------------------------

    def set_default_options(self, gamma=1, mu_s=1, alpha=0.1):
        """

        Default option for the integrator

        * The value of gamma (gyromagnetic ratio) for a free electron
        is 1.76e11

        """
        self.default_c = -1
        self._alpha[:] = alpha

        # When we create the simulation, mu_s is set to the default value. This
        # is overriden when calling the set_mu_s method from the Siulation
        # class or when setting mu_s directly (property)
        self._mu_s[:] = mu_s
        self.mu_s_const = mu_s
        self.gamma = gamma
        self.do_precession = True

    def set_tols(self, rtol=1e-8, atol=1e-10):
        """
        Set the relative and absolute tolerances for the CVODE integrator
        """
        self.integrator.set_options(rtol, atol)

    # Don't know the uselfulness of this function:
    # def set_options(self, rtol=1e-8, atol=1e-10):
    #     self.set_tols(rtol, atol)

    def compute_effective_field(self, t):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the micromagnetic Energy class
        """
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
        dm = (m1 - m0).reshape((-1, 3))
        max_dm = np.max(np.sqrt(np.sum(dm**2, axis=1)))
        max_dmdt = max_dm / dt
        return max_dmdt

    def stat(self):
        return self.integrator.stat()

    def reset_integrator(self, t=0):
        self.integrator.reset(self.spin, t)
        self.t = t  # also reinitialise the simulation time and step
        self.step = 0

    def sundials_rhs(self, t, y, ydot):
        """
        Defined in the corresponding driver class
        """
        pass

    def sundials_jtn(self, mp, Jmp, t, m, fy):
        # we can not copy mp to self.spin since m and self.spin is one object.
        #self.spin[:] = mp[:]
        print('NO jac...........')
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

    def run_until(self, t):

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
            return

        ode = self.integrator

        self.spin_last[:] = self.spin[:]

        flag = ode.run_until(t)

        if flag < 0:
            raise Exception("Run cython run_until failed!!!")

        self.spin[:] = ode.y[:]

        self.t = t
        self.step += 1

        # update field before saving data
        self.compute_effective_field(t)
        self.data_saver.save()

    def relax(self, dt=1e-11, stopping_dmdt=0.01,
              max_steps=1000, save_m_steps=100, save_vtk_steps=100):

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

            print('step=%d, time=%0.3g, max_dmdt=%0.3g ode_step=%0.3g'
                  % (self.step, self.t, dmdt, cvode_dt))

            if dmdt < stopping_dmdt:
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
        Save a VTK file with the magnetisation vector field and magnetic
        moments as cell data. Magnetic moments are saved in units of
        Bohr magnetons

        NOTE: It is recommended to use a *cell to point data* filter in
        Paraview or Mayavi to plot the vector field
        """
        self.vtk.save_vtk(self.spin.reshape(-1, 3),
                          self._mu_s / const.mu_B,
                          step=self.step
                          )

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
