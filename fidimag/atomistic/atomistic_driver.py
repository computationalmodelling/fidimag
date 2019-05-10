from __future__ import division
from __future__ import print_function

from fidimag.common.driver_base import DriverBase

import fidimag.extensions.a_clib as clib
import fidimag.extensions.cvode as cvode
from fidimag.common.vtk import VTK
import fidimag.common.constant as const

import numpy as np


class AtomisticDriver(DriverBase):
    """

    A class with shared methods and properties for different drivers to solve
    the Landau-Lifshitz-Gilbert equation

    Variables that are proper of the driver class:

        * alpha (damping)
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

        From DriverBase: t, spin_last, dm_dt, integrator_tolerances_set, step

    """

    def __init__(self, mesh, spin, mu_s, mu_s_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac,
                 integrator='sundials'
                 ):

        super(AtomisticDriver, self).__init__()

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
        self._pins = pins
        self.interactions = interactions
        # Strings are not referenced, this is a copy:
        self.name = name

        # The following are proper of the driver class: (see DriverBase) ------
        # See also the set_default_options() function

        self.n = self.mesh.n
        self.n_nonzero = self.mesh.n  # number of spins that are not zero
                                      # We check this in the set_Ms function

        self.initiate_variables(self.n)
        self.set_default_options()

        # Integrator options --------------------------------------------------

        # In the old code, it seemed that self.sundials_jtn was not defined
        # anywhere else. Here, we will use the same integrators than in the
        # micromagnetic code, where we have self.sundials_jtimes instead of jtn
        self.set_integrator(integrator, use_jac)

        self.set_tols()

        # When initialising the integrator in the self.integrator call, the
        # CVOde class calls the set_initial_value function (with flag_m=0),
        # which initialises a new integrator and allocates memory in this
        # process.  Now, when we set the magnetisation, we will use the same
        # memory setting this flag_m to 1, so instead of calling CVodeInit we
        # call CVodeReInit. If don't, memory is allocated in every call of
        # set_m
        self.flag_m = 1

        # Factor for the dmdt magnitude in the relaxation function
        self._dmdt_factor = 1.

        # Savers --------------------------------------------------------------

        # VTK saver for the magnetisation/spin field
        self.VTK = VTK(self.mesh,
                       directory='{}_vtks'.format(self.name),
                       filename='m'
                       )

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
        # is overriden when calling the set_mu_s method from the Simulation
        # class or when setting mu_s directly (property)
        # David Tue 19 Jun 2018: Not a very clear thing to do, we must set a
        # WARNING
        self._mu_s[:] = mu_s
        self._mu_s_inv[:] = 1. / mu_s

        self.mu_s_const = mu_s
        self.gamma = gamma
        self.do_precession = True

    # Don't know the uselfulness of this function:
    # def set_options(self, rtol=1e-8, atol=1e-10):
    #     self.set_tols(rtol, atol)

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
        self.VTK.reset_data()

        # Here we save both Ms and spins as cell data
        self.VTK.save_scalar(self._mu_s / const.mu_B, name='mu_s')
        self.VTK.save_vector(self.spin.reshape(-1, 3), name='spins')

        self.VTK.write_file(step=self.step)
