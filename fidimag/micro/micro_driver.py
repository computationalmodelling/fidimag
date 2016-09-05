from __future__ import division
from __future__ import print_function

from fidimag.common.driver_base import DriverBase

import numpy as np
from fidimag.common.integrators import CvodeSolver, CvodeSolver_OpenMP, StepIntegrator
# from fidimag.common.fileio import DataSaver, DataReader
from fidimag.common.vtk import VTK
import time


class MicroDriver(DriverBase):
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

        super(MicroDriver, self).__init__()

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

        # Savers --------------------------------------------------------------

        # VTK saver for the magnetisation/spin field
        self.VTK = VTK(self.mesh,
                       directory='{}_vtks'.format(self.name),
                       filename='m'
                       )

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

        # OOMMF convention is to check if the spins have moved by ~1 degree in
        # a nanosecond in order to stop a simulation, so we set this scale for
        # dm/dt

        # ONE_DEGREE_PER_NANOSECOND:
        self._dmdt_factor = (2 * np.pi / 360) / 1e-9

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

    # -------------------------------------------------------------------------
    # Save functions ----------------------------------------------------------
    # -------------------------------------------------------------------------

    def save_vtk(self):
        """
        Save a VTK file with the magnetisation vector field (vector data)
        and the saturation magnetisation values (scalar data) as
        cell data
        """
        self.VTK.reset_data()

        # Here we save both Ms and spins as cell data
        self.VTK.save_scalar(self._Ms, name='mu_s')
        self.VTK.save_vector(self.spin.reshape(-1, 3), name='spins')

        self.VTK.write_file(step=self.step)
