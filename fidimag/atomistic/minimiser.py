from __future__ import division
import numpy as np
import fidimag.extensions.clib as clib
import fidimag.common.helper as helper
import fidimag.common.constant as const

from .atomistic_driver import AtomisticDriver


class Minimiser(AtomisticDriver):
    """

    A simple minimisation algorithm

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
        super(Minimiser, self).__init__(mesh, spin, mu_s, mu_s_inv, field,
                                        pins, interactions, name,
                                        data_saver,
                                        use_jac=use_jac,
                                        integrator=integrator
                                        )

        self.t = 1e-4
        self.alpha = 0.1
        self.spin_last = np.zeros_like(spin)
        self._new_spin = np.zeros_like(spin)

        # self.set_options()

    def run_step(self):

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()

        self._new_spin[self._material] = (self.spin + self.field)[self._material]
        clib.normalise_spin(self._new_spin, self._pins, self.n)

        self._new_spin[self._material] = (self.spin_last +
                                          self.alpha_field * (self._new_spin - self.spin_last))[self._material]
        self.spin[:] = self._new_spin[:]
        clib.normalise_spin(self.spin, self._pins, self.n)

    def update_effective_field(self):

        self.field[:] = 0
        for obj in self.interactions:
            self.field += obj.compute_field(t=0, spin=self.spin)

    def minimise(self, stopping_dm=1e-2, max_steps=2000, save_data_steps=10):
        """

        """

        self.step = 0

        self.alpha_field = np.repeat(self.alpha, 3)

        # Only update site with mu_s > 0 which are not pinned
        self._material = np.logical_and(np.repeat(self._mu_s / const.mu_B, 3) > 1e-10,
                                        np.repeat(1 - self._pins, 3).astype(np.bool))

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()
        while self.step < max_steps:

            self.run_step()

            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis=1)))
            print("#max_dm={:<10.3g} step={}".format(max_dm, self.step))

            if max_dm < stopping_dm and self.step > 0:
                break

            if self.step % save_data_steps == 0:
                # update field before saving data
                self.update_effective_field()
                self.data_saver.save()

            self.step += 1

    def relax(self):
        print('Not implemented for the minimizer')
