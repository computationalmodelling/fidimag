from __future__ import division
import numpy as np
import fidimag.extensions.clib as clib
# import fidimag.common.constant as const

from .minimiser_base import MinimiserBase


class SimpleMinimiser(MinimiserBase):
    """
    A simple minimisation algorithm, where the evolution of the magnetisation
    follows the system's torque

        CHECK:

        h = (m_i + H) / || (m_i + H) ||
        m_i+1 = m_i + alpha * (h_i - m_i)
    """

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        # Inherit from the base minimiser class
        super(SimpleMinimiser, self).__init__(mesh, spin,
                                              magnetisation, magnetisation_inv,
                                              field,
                                              pins,
                                              interactions,
                                              name,
                                              data_saver
                                              )

        self.t = 1e-4
        self._alpha = 0.1
        self._alpha_field = self._alpha * np.ones_like(self.spin)
        self.spin_last = np.zeros_like(spin)
        self._new_spin = np.zeros_like(spin)

    @property
    def alpha(self):
        """
        Returns the array with the spatially dependent Gilbert damping
        per mesh/lattice site
        """
        return self._alpha

    @alpha.setter
    def set_alpha(self, value):
        """
        """
        self._alpha = value
        self._alpha_field =  value * np.ones_like(self.spin)

    def run_step(self):

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()

        self._new_spin[self._material] = (self.spin + self.field)[self._material]
        clib.normalise_spin(self._new_spin, self._pins, self.n)

        self._new_spin[self._material] = (self.spin_last +
                                          self._alpha_field * (self._new_spin - self.spin_last))[self._material]
        self.spin[:] = self._new_spin[:]
        clib.normalise_spin(self.spin, self._pins, self.n)

    def minimise(self, stopping_dm=1e-2, max_steps=2000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1000
                 ):
        """

        """

        self.step = 0

        # Only update site with magnetisation > 0 which are not pinned
        self._material = np.logical_and(np.repeat(self._magnetisation, 3) > 0.0,
                                        np.repeat(1 - self._pins, 3).astype(np.bool))

        self.spin_last[:] = self.spin[:]
        self.compute_effective_field()
        while self.step < max_steps:

            self.run_step()

            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis=1)))

            if self.step % log_steps == 0:
                print("#max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
                    np.max(np.abs(self.tau)),
                    max_dm, self.step))

            if max_dm < stopping_dm and self.step > 0:
                print("FINISHED AT: max_tau={:<8.3g} max_dm={:<10.3g} counter={}".format(
                      np.max(np.abs(self.tau)),
                      max_dm, self.step))

                self.compute_effective_field()
                self.data_saver.save()

                break

            if self.step % save_data_steps == 0:
                # update field before saving data
                self.compute_effective_field()
                self.data_saver.save()

            if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                self.save_vtk()
            if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                self.save_m()

            self.step += 1
