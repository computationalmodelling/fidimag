from __future__ import division
import numpy as np
import fidimag.extensions.common_clib as clib
# Change int he future to common clib:
import fidimag.extensions.clib as atom_clib
import sys
from .minimiser_base import MinimiserBase


class SteepestDescent(MinimiserBase):
    """

    This class is the driver to minimise a system using an optimised Steepest
    Descent algorithm defined in [1].

    The evolution step is written as

        m_i+1 = FM * (FP * m_i - 4 * tau * (m_i X m_i x H))

    where

        FM = (1 - tau^2 * (m_i X H)^2)
        FP = (1 + tau^2 * (m_i X H)^2)

    and tau is a "time step" that needs to be defined according to
    eq. 10 of [1].

    NOTES:

    - We use the simplest criteria for choosing tau. However, it is not exactly
      clear how to deal with the denominators defined in tau when they are
      zero.

      We are not using these criteria but might be necessary:
      For now we refer to the methods defined in the minimizer branch of
      https://github.com/MicroMagnum/MicroMagnum.

    - The effective field H is defined in Tesla units (it works in the
      atomistic case). This is defined in the MuMax3 code, so we scale the
      field with mu0 when using this class with the micromagnetic classes.

    - Methods for the calculations using Numpy are defined for debugging. In
      this class we use a C library to speed up the evolution step.


    REFS:

    [1] Exl et al., Journal of Applied Physics 115, 17D118 (2014).
    https://doi.org/10.1063/1.4862839


    """

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        # Define
        super(SteepestDescent, self).__init__(mesh, spin,
                                              magnetisation, magnetisation_inv,
                                              field,
                                              pins,
                                              interactions,
                                              name,
                                              data_saver)

        # ---------------------------------------------------------------------
        # Variables defined in this SteepestDescent

        self.mxH = np.zeros_like(self.field)
        self.mxmxH = np.zeros_like(self.field)
        self.mxmxH_last = np.zeros_like(self.field)
        # time as zero
        self.t = 0
        self.tau = 1e-4
        # self._n_field = np.zeros_like(self.field)

        # Threshold values for the criteria to choose the T factor.
        # They are defined as properties with the corr decoration
        self._tmax = 1e-1
        self._tmin = 1e-16

        # Scaling of the field
        self.scale = 1.

    @property
    def tmax(self):
        return self._tmax

    @tmax.setter
    def tmax(self, t):
        self._tmax = t
        # self._tmax_arr = t * np.ones((len(self.tau), 2))

    @property
    def tmin(self):
        return self._tmin

    @tmin.setter
    def tmin(self, t):
        self._tmin = t
        # self._tmin_arr = t * np.ones((len(self.tau), 2))

    # Same as:
    # tmin = property(fget=set_tmin, fset=set_tmin)

    def normalise_field(self, a):
        norm = np.sqrt(np.sum(a.reshape(-1, 3) ** 2, axis=1))
        norm_a = a.reshape(-1, 3) / norm[:, np.newaxis]
        norm_a.shape = (-1,)
        return norm_a

    def field_cross_product(self, a, b):
        aXb = np.cross(a.reshape(-1, 3), b.reshape(-1, 3))
        return aXb.reshape(-1,)

    def run_step(self):
        """
        Numpy version of the calculation
        """

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
                    - 4 * self.tau * self.mxmxH
                    # this term should be zero:
                    # + (2 * (self.tau ** 2) * np.sum(self.mxH * self.spin, axis=1))[:, np.newaxis] * self.mxH
                    )
        new_spin = new_spin / factor_plus[:, np.newaxis]

        self.mxH.shape = (-1,)
        self.mxmxH.shape = (-1,)
        self.spin.shape = (-1,)
        new_spin.shape = (-1,)

        self.spin_last[:] = self.spin[:]
        self.spin[:] = new_spin[:]
        atom_clib.normalise_spin(self.spin, self._pins, self.n)

        # ---------------------------------------------------------------------

        # Update the effective field, torques and time step for the next iter
        self.compute_effective_field()

        self.mxmxH_last[:] = self.mxmxH[:]
        self.mxH[:] = self.field_cross_product(self.spin, self.scale * self.field)[:]
        self.mxmxH[:] = self.field_cross_product(self.spin, self.mxH)[:]

        # ---------------------------------------------------------------------
        # Define the time step tau

        ds = (self.spin - self.spin_last).reshape(-1, 3)
        dy = (self.mxmxH - self.mxmxH_last).reshape(-1, 3)

        if self.step % 2 == 0:
            num = np.sum(ds * ds)
            den = np.sum(ds * dy)
        else:
            num = np.sum(ds * dy)
            den = np.sum(dy * dy)

        # The criteria is taken from the Micromagnum code
        if den == 0:
            self.tau = self._tmax
        else:
            self.tau = num / den

        # Set the minimum between the abs value of tau and the max tolerance
        # self.tau = np.sign(self.tau) * np.max(self._tmin_arr, axis=1)
        # Set the maximum between the previous minimum and the min tolerance
        # self.tau = np.sign(self.tau) * max(min(np.abs(self.tau),
        #                                        self._tmax),
        #                                    self._tmin)

        # ---------------------------------------------------------------------

    def run_step_CLIB(self):
        """
        C version of the calculation from common/lib/steepest_descent.c
        """

        clib.compute_sd_spin(self.spin, self.spin_last,
                             self._magnetisation,
                             self.mxH, self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n
                             )

        self.compute_effective_field()

        # Notice that the field is scaled (in the micro class we use Tesla)
        clib.compute_sd_step(self.spin, self.spin_last,
                             self._magnetisation,
                             self.scale * self.field,
                             self.mxH, self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n, self.step,
                             self._tmin, self._tmax
                             )

    def minimise(self, stopping_dm=1e-3, max_steps=5000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_every=1000, printing=True,
                 initial_t_step=1e-2
                 ):
        """
        Run the minimisation until meeting the stopping_dm criteria
        """

        # Rewrite tmax and tmin arrays and variable
        self.tmax = self._tmax
        self.tmin = self._tmin

        self.step = 0
        # Initial "time" step: the algorithm seems sensitive to this value
        self.tau = initial_t_step

        self.spin_last[:] = self.spin[:]
        self.compute_effective_field()
        self.mxH[:] = self.field_cross_product(self.spin, self.scale * self.field)[:]
        self.mxmxH[:] = self.field_cross_product(self.spin, self.mxH)[:]
        self.mxmxH_last[:] = self.mxmxH[:]
        while self.step < max_steps:
            self.run_step_CLIB()
            # Vectorised calculation with Numpy:
            # self.run_step()
            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis=1)))

            if printing:
                if self.step % log_every == 0:
                    # print("#{:<4} t={:<8.3g} dt={:.3g} max_dmdt={:.3g}
                    print("#{:<4} max_tau={:<8.3g} max_dm={:<10.3g}".format(self.step,
                            np.max(np.abs(self.tau)),
                            max_dm))

            if max_dm < stopping_dm and self.step > 0:
                print("#{:<4} max_tau={:<8.3g} max_dm={:<10.3g}".format(self.step,
                                                                        np.max(np.abs(self.tau)),
                                                                        max_dm)
                      )
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

        if self.step == max_steps:
            sys.stderr.write("Warning: minimise did not converge in {} steps - maxdm = {}".format(self.step, max_dm))
