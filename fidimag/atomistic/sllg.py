from __future__ import division
import numpy as np
import fidimag.extensions.clib as clib
import fidimag.common.helper as helper
from .llg import LLG
import fidimag.common.constant as const


class SLLG(LLG):

    def __init__(self, mesh, name='unnamed'):
        """Simulation object.

        *Arguments*

          name : the Simulation name (used for writing data files, for examples)

        """
        super(SLLG, self).__init__(mesh, name=name)

        self._T = np.zeros(self.n, dtype=np.float)
        self.next_spin = np.zeros(3*self.n, dtype=np.float)
        self.eta = np.zeros(3*self.n, dtype=np.float)
        self.dm1 = np.zeros(3*self.n, dtype=np.float)
        self.dm2 = np.zeros(3*self.n, dtype=np.float)

        self.set_options()

    def get_T(self):
        return self._T

    def set_T(self, T0):
        self._T[:] = helper.init_scalar(T0, self.mesh)

    T = property(get_T, set_T)

    def set_options(self, dt=1e-15, theta=1.0, gamma=const.gamma, k_B=const.k_B, seed=100):

        clib.init_random(seed)
        self.gamma = gamma
        self.k_B = k_B
        self.dt = dt
        self.theta = theta
        self.theta1 = 1-0.5/theta;
	self.theta2 = 0.5/theta;

    def run_step(self):
        
        clib.random_number_array(self.eta)
        
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

        self.step += 1
        self.t = self.dt*self.step

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

        clib.normalise_spin(self.spin, self.n)
        


    def update_effective_field(self, y, t):

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t, spin=y)


    def run_until(self, t):

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.saver.save()
            return

        self.spin_last[:] = self.spin[:]

        while (self.t<t):
            self.run_step()

        # update field before saving data
        self.compute_effective_field(t)
        self.saver.save()
