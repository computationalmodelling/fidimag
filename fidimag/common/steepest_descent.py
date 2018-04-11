from __future__ import division
import numpy as np
import fidimag.extensions.common_clib as clib
import fidimag.common.helper as helper
import fidimag.common.constant as const
from fidimag.common.vtk import VTK

from .driver_base import DriverBase


class SteepestDescent(object):
    """

    This class is the driver to minimise a system using a Steepest Descent
    algorithm

    NOTE: We are inheriting from DriverBase, but it would be better if we
          create a MinimiserBase class in the future, removing the
          methods from the Driver class that are not being used at all

    """

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 use_jac=False,
                 integrator=None
                 ):

        # ---------------------------------------------------------------------
        # These are (ideally) references to arrays taken from the Simulation
        # class. Variables with underscore are arrays changed by a property in
        # the simulation class
        self.mesh = mesh
        self.spin = spin

        # A reference to either mu_s or Ms to use a common lib for this
        # minimiser
        self._magnetisation = magnetisation
        self._magnetisation_inv = magnetisation_inv

        self.field = field
        self._pins = pins
        self.interactions = interactions
        # Strings are not referenced, this is a copy:
        self.name = name

        self.data_saver = data_saver

        # ---------------------------------------------------------------------
        # Variables defined in this class

        self.spin_last = np.ones_like(self.spin)
        self.n = self.mesh.n

        # VTK saver for the magnetisation/spin field
        self.VTK = VTK(self.mesh,
                       directory='{}_vtks'.format(self.name),
                       filename='m'
                       )

        self.mxH = np.zeros_like(self.field)
        self.mxmxH = np.zeros_like(self.field)
        self.mxmxH_last = np.zeros_like(self.field)
        self.t = 1e-4
        self.tau = 1e-4 * np.ones(len(self.spin) // 3)
        # self._n_field = np.zeros_like(self.field)

        # self.set_options()
        self._tmax = 1e-2
        self._tmin = 1e-16

        # If the magnetisation changes, this needs updating! 
        self.scale = np.ones_like(self._magnetisation)

    @property
    def tmax(self):
        return self._tmax

    @tmax.setter
    def tmax(self, t):
        self._tmax = t
        self.__tmax = t * np.ones((len(self.tau), 2))

    @property
    def tmin(self):
        return self._tmin

    @tmin.setter
    def tmin(self, t):
        self._tmin = t
        self.__tmin = t * np.ones((len(self.tau), 2))

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

        self.spin_last[:] = self.spin[:]
        self.spin[:] = new_spin[:]
        clib.normalise_spin(self.spin, self._pins, self.n)

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
        self.tau = self._tmax * np.ones_like(num)
        self.tau[den != 0] = num[den != 0] / den[den != 0]

        tau_signs = np.sign(self.tau)

        self.__tmax[:, 0] = np.abs(self.tau)
        # Set the minimum between the abs value of tau and the max tolerance
        self.__tmin[:, 0] = np.min(self.__tmax, axis=1)
        # Set the maximum between the previous minimum and the min tolerance
        self.tau = tau_signs * np.max(self.__tmin, axis=1)

        # ---------------------------------------------------------------------

        # clib.normalise_spin(self.spin, self._pins, self.n)
        # self.spin[:] = self.normalise_field(self.spin)[:]
        # self.compute_rhs(self.tau)

    def run_step_CLIB(self):
        """
        Only use when run_step with Numpy is working
        """

        clib.compute_sd_spin(self.spin, self.spin_last,
                             self.mxH, self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n
                             )

        self.update_effective_field()

        clib.compute_sd_step(self.spin, self.spin_last,
                             self.field, self.scale,
                             self.mxH, self.mxmxH, self.mxmxH_last,
                             self.tau, self._pins,
                             self.n, self.step,
                             self._tmin, self._tmax
                             )

    def update_effective_field(self):

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t=0, spin=self.spin)

    def minimise(self, stopping_dm=1e-2, max_steps=2000, 
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1000
                 ):


        # Rewrite tmax and tmin arrays and variable
        self.tmax = self._tmax
        self.tmin = self._tmin

        self.step = 0

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()
        # self._n_field[:] = self.field[:]
        # clib.normalise_spin(self._n_field, self._pins, self.n)
        self.mxH[:] = self.field_cross_product(self.spin, self.field)[:]
        self.mxmxH[:] = self.field_cross_product(self.spin, self.mxH)[:]
        self.mxmxH_last[:] = self.mxmxH[:]
        while self.step < max_steps:

            self.run_step_CLIB()

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

                self.update_effective_field()
                self.data_saver.save()

                break

            # print('spin=', self.spin.reshape(-1, 3)[33])
            # print('field=', self.field.reshape(-1, 3)[33])
            # print('tau', self.tau[33])

            if self.step % save_data_steps == 0:
                # update field before saving data
                self.update_effective_field()
                self.data_saver.save()

            if (save_vtk_steps is not None) and (self.step % save_vtk_steps == 0):
                self.save_vtk()
            if (save_m_steps is not None) and (self.step % save_m_steps == 0):
                self.save_m()

            self.step += 1

    def relax(self):
        print('Not implemented for the SD minimiser')

    # -------------------------------------------------------------------------

    def compute_effective_field(self, t):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the corresponding Energy class
        """

        # self.spin[:] = y[:]

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t)

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
        self.VTK.save_scalar(self._magnetisation / const.mu_B, name='magnetisation')
        self.VTK.save_vector(self.spin.reshape(-1, 3), name='spins')

        self.VTK.write_file(step=self.step)

    def save_m(self, ZIP=False):
        """
        Save the magnetisation/spin vector field as a numpy array in
        a NPY file. The files are saved in the `{name}_npys` folder, where
        `{name}` is the simulation name, with the file name `m_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """

        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)
        if ZIP:
            with zipfile.ZipFile('%s_m.zip'%self.name, 'a') as myzip:
                myzip.write(name)
            try:
                os.remove(name)
            except OSError:
                pass

    def save_skx(self):
        """
        Save the skyrmion number density (sk number per mesh site)
        as a numpy array in a NPY file.
        The files are saved in the `{name}_skx_npys` folder, where
        `{name}` is the simulation name, with the file name `skx_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """
        if not os.path.exists('%s_skx_npys' % self.name):
            os.makedirs('%s_skx_npys' % self.name)
        name = '%s_skx_npys/m_%g.npy' % (self.name, self.step)

        # The _skx_number array is defined in the SimBase class in Common
        np.save(name, self._skx_number)
