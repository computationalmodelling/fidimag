from __future__ import division
import numpy as np
import fidimag.extensions.clib as clib
import fidimag.common.helper as helper
import fidimag.common.constant as const

from .atomistic_driver import AtomisticDriver


class SteepestDescent(AtomisticDriver):
    """

    This class is the driver to minimise a system using a Steepest Descent
    algorithm



    This class inherits common methods to evolve the system using CVODE, from
    the micro_driver.AtomisticDriver class. Arrays with the system information
    are taken as references from the main micromagnetic Simulation class

    """

    def __init__(self, mesh, spin, mu_s, mu_s_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        # self.mesh = mesh
        # self.spin = spin
        # self.mu_s mu_s
        # self._mu_s_inv = mu_s_inv
        # self.field = field
        # self.pins = pins
        # self.interactions = interactions
        # self.name = name
        # self.data_saver = data_saver

        # Inherit from the driver class
        super(SteepestDescent, self).__init__(mesh, spin, mu_s, mu_s_inv, field,
                                              pins, interactions, name,
                                              data_saver,
                                              use_jac=use_jac,
                                              integrator=integrator
                                              )

        self.mxH = np.zeros_like(self.field)
        self.mxmxH = np.zeros_like(self.field)
        self.mxmxH_last = np.zeros_like(self.field)
        self.t = 1e-4
        self.tau = 1e-4 * np.ones(len(self.spin) // 3)
        # self._n_field = np.zeros_like(self.field)

        # self.set_options()

    def normalise_field(self, a):
        norm = np.sqrt(np.sum(a.reshape(-1, 3) ** 2, axis=1))
        norm_a = a.reshape(-1, 3) / norm[:, np.newaxis]
        norm_a.shape = (-1,)
        return norm_a

    def field_cross_product(self, a, b):
        aXb = np.cross(a.reshape(-1, 3), b.reshape(-1, 3))
        return aXb.reshape(-1,)

    def run_step(self):

        # ---------------------------------------------------------------------

        self.mxH.shape = (-1, 3)
        self.mxmxH.shape = (-1, 3)
        self.spin.shape = (-1, 3)

        mxH_sq_norm = np.sum(self.mxH ** 2, axis=1)
        factor_plus = 4 + (self.tau ** 2) * mxH_sq_norm
        factor_minus = 4 - (self.tau ** 2) * mxH_sq_norm

        # Compute: m[i+1] = ((4 - t^2 A^2) * m[i] - 4 * t * m[i] x m[i] x H) / (4 + t^2 A^2) 
        # where "t = self.tau" is the time step and "A = m[i] x H"
        new_spin = (factor_minus[:, np.newaxis] * self.spin
                    - (4 * self.tau)[:, np.newaxis] * self.mxmxH
                    # this term should be zero:
                    # + (2 * (self.tau ** 2) * np.sum(self.mxH * self.spin, axis=1))[:, np.newaxis] * self.mxH
                    )
        new_spin = new_spin / factor_plus[:, np.newaxis]

        self.mxH.shape = (-1,)
        self.mxmxH.shape = (-1,)
        self.spin.shape = (-1,)
        new_spin.shape = (-1,)

        clib.normalise_spin(self.spin, self._pins, self.n)
        self.spin_last[:] = self.spin[:]
        self.spin[:] = new_spin[:]

        # ---------------------------------------------------------------------

        # Update the effective field, torques and time step for the next iter
        self.update_effective_field()
        # self._n_field[:] = self.field[:]
        # clib.normalise_spin(self._n_field, self._pins, self.n)

        self.mxmxH_last[:] = self.mxmxH[:]
        self.mxH[:] = self.field_cross_product(self.spin, self.field)[:]
        self.mxmxH[:] = self.field_cross_product(self.spin, self.mxH)[:]

        # clib.normalise_spin(self.spin, self._pins, self.n)

        # ---------------------------------------------------------------------
        # self.tau = self.get_time_step()
        ds = (self.spin - self.spin_last).reshape(-1, 3)
        dy = (self.mxmxH - self.mxmxH_last).reshape(-1, 3)

        if self.step % 2 == 0:
            num = np.sum(ds * ds, axis=1)
            den = np.sum(ds * dy, axis=1)
        else:
            num = np.sum(ds * dy, axis=1)
            den = np.sum(dy * dy, axis=1)
      
        # Set to tmin
        self.tau = 1e-2 * np.ones_like(num)
        self.tau[den != 0] = num[den != 0] / den[den != 0]

        tau_signs = np.sign(self.tau)

        self._tmax[:, 0] = np.abs(self.tau)
        # Set the minimum between the abs value of tau and the max tolerance
        self._tmin[:, 0] = np.min(self._tmax, axis=1)
        # Set the maximum between the previous minimum and the min tolerance
        self.tau = tau_signs * np.max(self._tmin, axis=1)

        # ---------------------------------------------------------------------

        # clib.normalise_spin(self.spin, self._pins, self.n)
        # self.spin[:] = self.normalise_field(self.spin)[:]
        # self.compute_rhs(self.tau)

    def run_step_CLIB(self):
        """
        Only use when run_step with Numpy is working
        """

        clib.compute_sd_spin(self.spin, self.spin_last,
                             self.field, self.mxH,
                             self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n)

        self.update_effective_field()

        clib.compute_sd_step(self.spin, self.spin_last,
                             self.field, self.mxH,
                             self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n, self.step)

    def update_effective_field(self):

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t=0, spin=self.spin)

    def minimise(self, stopping_dm=1e-2, max_steps=2000, 
                 tmax=1e-2, tmin=1e-16):
        self.step = 0

        self._tmax = tmax * np.ones((len(self.tau), 2))
        self._tmin = tmin * np.ones((len(self.tau), 2))

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()
        # self._n_field[:] = self.field[:]
        # clib.normalise_spin(self._n_field, self._pins, self.n)
        self.mxH[:] = self.field_cross_product(self.spin, self.field)[:]
        self.mxmxH[:] = self.field_cross_product(self.spin, self.mxH)[:]
        self.mxmxH_last[:] = self.mxmxH[:]
        while self.step < max_steps:

            self.run_step()

            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis=1)))
            print("#max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
                np.max(np.abs(self.tau)),
                max_dm, self.step))
            if max_dm < stopping_dm and self.step > 0:
                break

            self.step += 1

            # print('spin=', self.spin.reshape(-1, 3)[33])
            # print('field=', self.field.reshape(-1, 3)[33])
            # print('tau', self.tau[33])

            # update field before saving data
            # self.update_effective_field()
            # self.data_saver.save()

        # clib.normalise_spin(self.spin, self._pins, self.n)

    def relax(self):
        print('Not implemented for the SD minimiser')
