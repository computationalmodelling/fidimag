from __future__ import division
import numpy as np
import fidimag.extensions.clib as clib
import fidimag.common.helper as helper
import fidimag.common.constant as const

from .atomistic_driver import AtomisticDriver


class SLLG(AtomisticDriver):
    """

    This class is the driver to solve the Stochastic Landau Lifshitz Gilbert
    equation.  The equation is given by:


          ds        -gamma
         ---- =    --------  ( s X H_eff  ...
          dt             2
                  ( 1 + a  )


    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.AtomisticDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

    """

    def __init__(self, mesh, spin, mu_s, mu_s_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator='sundials'
                 ):

        # Inherit from the driver class
        super(SLLG, self).__init__(mesh, spin, mu_s, mu_s_inv, field,
                                   pins, interactions, name,
                                   data_saver,
                                   use_jac=use_jac,
                                   integrator=integrator
                                   )

        self._T = np.zeros(self.n, dtype=np.float64)
        self.next_spin = np.zeros(3*self.n, dtype=np.float64)
        self.eta = np.zeros(3*self.n, dtype=np.float64)
        self.dm1 = np.zeros(3*self.n, dtype=np.float64)
        self.dm2 = np.zeros(3*self.n, dtype=np.float64)

        self.minor_step = 0
        self.mt19937 = clib.rng_mt19937()

        self.set_options()

    def get_T(self):
        return self._T

    def set_T(self, T0):
        self._T[:] = helper.init_scalar(T0, self.mesh)

    T = property(get_T, set_T)

    def set_options(self, dt=1e-15, theta=1.0,
                    gamma=const.gamma,
                    k_B=const.k_B,
                    seed=100
                    ):

        self.mt19937.set_seed(seed)
        self.gamma = gamma
        self.k_B = k_B
        self.dt = dt
        self.theta = theta
        self.theta1 = 1-0.5/theta
        self.theta2 = 0.5/theta

    def run_step(self):

        self.mt19937.fill_vector_gaussian(self.eta)

        #step1
        self.update_effective_field(self.spin, self.t)
        clib.compute_llg_rhs_dw(self.dm1,
                                self.spin,
                                self.field,
                                self._T,
                                self._alpha,
                                self._mu_s_inv,
                                self.eta,
                                self._pins,
                                self.n,
                                self.gamma,
                                self.dt)
        self.next_spin = self.spin + self.theta * self.dm1

        self.minor_step += 1
        self.t = self.dt*self.minor_step

        #step2
        self.update_effective_field(self.next_spin, self.t)
        clib.compute_llg_rhs_dw(self.dm2,
                                self.next_spin,
                                self.field,
                                self._T,
                                self._alpha,
                                self._mu_s_inv,
                                self.eta,
                                self._pins,
                                self.n,
                                self.gamma,
                                self.dt)
        self.spin += (self.theta1*self.dm1+self.theta2*self.dm2)

        clib.normalise_spin(self.spin, self._pins, self.n)

    def update_effective_field(self, y, t):

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t, spin=y)

    def run_until(self, t):

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
            return

        self.spin_last[:] = self.spin[:]

        while (self.t < t):
            self.run_step()
        self.step += 1

        # update field before saving data
        self.compute_effective_field(t)
        self.data_saver.save()
