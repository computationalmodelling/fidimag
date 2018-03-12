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

        # self.set_options()

    def field_cross_product(self, a, b):
        aXb = np.cross(a.reshape(-1, 3), b.reshape(-1, 3))
        return aXb.reshape(-1,)

    def get_time_step(self):
        ds = self.spin - self.spin_last
        dy = self.mxmxH - self.mxmxH_last

        if self.counter % 2 == 0:
            num = np.dot(ds, ds)
            den = np.dot(ds, dy)
        else:
            num = np.dot(ds, dy)
            den = np.dot(dy, dy)

        if den != 0:
            return num / den
        else:
            return 1e-4

    def compute_rhs(self):

        mxH_sq_norm = np.sum(self.mxH ** 2)

        factor_plus = 4 + (self.t ** 2) * mxH_sq_norm
        factor_minus = 4 - (self.t ** 2) * mxH_sq_norm

        self.spin = factor_minus * self.spin - 4 * self.t * self.mxmxH

        clib.normalise_spin(self.spin, self._pins, self.n)

        self.spin /= factor_plus

    def run_step(self):

        self.update_effective_field()
        self.mxH = self.field_cross_product(self.spin, self.field)
        self.mxmxH = self.field_cross_product(self.spin, self.mxH)
        self.t = self.get_time_step()
        self.compute_rhs()

    def update_effective_field(self):

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t=0, spin=self.spin)

    def minimize(self, stopping_dmdt=1e-2, max_count=2000):
        self.counter = 0

        while self.counter < max_count:
            self.spin_last[:] = self.spin[:]
            self.mxmxH_last[:] = self.mxmxH[:]
            self.run_step()

            dmdt = self.compute_dmdt(self.t)
            print("#step={:<8.3g} max_dmdt={:<10.3g} counter={}".format(
                self.t,
                dmdt, self.counter))
            if dmdt < stopping_dmdt:
                break

            self.counter += 1

            # update field before saving data
            self.update_effective_field()
            self.data_saver.save()

    def relax(self):
        print('Not implemented for the minimizer')
