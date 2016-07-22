import numpy as np
from fidimag.common.integrators import SundialsIntegrator, StepIntegrator
from fidimag.common.fileio import DataSaver, DataReader
from fidimag.common.save_vtk import SaveVTK


class MicroDriver(object):
    """

    A class with shared properties for different drivers to solve the
    Landau-Lifshitz-Gilbert equation
    
    """

    def __init__(self, spin, alpha, field, pins,
                 integrator='sundials', 
                 use_jac=False
                 ):

        self.t = 0
        self._alpha = np.zeros(self.n, dtype=np.float)
        self._Ms_inv = np.zeros(self.n, dtype=np.float)
        self.spin_last = np.ones(3 * self.n, dtype=np.float)
        self.dm_dt = np.zeros(3 * self.n, dtype=np.float)

        self.interactions = []
        self.integrator_tolerances_set = False
        self.step = 0

        self.vtk = SaveVTK(self.mesh, name=name)
        self.set_default_options()

        if integrator == "sundials" and use_jac:
            self.integrator = SundialsIntegrator(self.spin,
                                                 self.sundials_rhs,
                                                 self.sundials_jtimes
                                                 )
        elif integrator == "sundials_diag":
            self.integrator = SundialsIntegrator(self.spin,
                                                 self.sundials_rhs,
                                                 linear_solver="diag"
                                                 )
        elif integrator == "sundials":
            self.integrator = SundialsIntegrator(self.spin,
                                                 self.sundials_rhs
                                                 )
        elif integrator == "euler" or integrator == "rk4":
            self.integrator = StepIntegrator(self.spin,
                                             self.step_rhs,
                                             integrator
                                             )
        else:
            raise NotImplemented("integrator must be sundials, euler or rk4")

        self.saver = DataSaver(self, name + '.txt')

        self.saver.entities['E_total'] = {
            'unit': '<J>',
            'get': lambda sim: sim.compute_energy(),
            'header': 'E_total'}

        self.saver.entities['m_error'] = {
            'unit': '<>',
            'get': lambda sim: sim.compute_spin_error(),
            'header': 'm_error'}

        self.saver.entities['skx_num'] = {
            'unit': '<>',
            'get': lambda sim: sim.skyrmion_number(),
            'header': 'skx_num'}

        self.saver.entities['rhs_evals'] = {
            'unit': '<>',
            'get': lambda sim: self.integrator.rhs_evals(),
            'header': 'rhs_evals'}

        self.saver.entities['real_time'] = {
            'unit': '<s>',
            'get': lambda _: time.time(),  # seconds since epoch
            'header': 'real_time'}

        self.saver.update_entity_order()

    def set_default_options(self, gamma=2.21e5, Ms=8.0e5, alpha=0.1):
        self.default_c = 1e11
        self._alpha[:] = alpha
        self._Ms[:] = Ms
        self.gamma = gamma
        self.do_precession = True

    def reset_integrator(self, t=0):
        self.integrator.reset(self.spin, t)
        self.t = t # also reinitialise the simulation time and step
        self.step = 0

    def set_tols(self, rtol=1e-8, atol=1e-10, max_ord=None, reset=True):
        if max_ord is not None:
            self.integrator.set_options(rtol=rtol, atol=atol, max_ord=max_ord)
        else:
            # not all integrators have max_ord (only VODE does)
            # and we don't want to encode a default value here either
            self.integrator.set_options(rtol=rtol, atol=atol)
        if reset:
            self.reset_integrator(self.t)

    def compute_effective_field(self, t):

        #self.spin[:] = y[:]

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t)

    def compute_effective_field_jac(self, t, spin):
        self.field[:] = 0
        for obj in self.interactions:
            if obj.jac:
                self.field += obj.compute_field(t, spin=spin)


