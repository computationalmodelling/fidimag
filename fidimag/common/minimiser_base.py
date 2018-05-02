from __future__ import division
import numpy as np
import os
import zipfile
import fidimag.common.constant as const
from fidimag.common.vtk import VTK


class MinimiserBase(object):
    """
    Base class for minimiser class. No dependency on CVODE
    """

    def __init__(self, mesh, spin,
                 magnetisation, magnetisation_inv, field, pins,
                 interactions,
                 name,
                 data_saver
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

        self.scale = 1.

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
        Python implementation of the step
        """
        pass

    def run_step_CLIB(self):
        """
        Cython implementation of the step. Normally you would define
        functions called in this method, in the lib/ folder
        """
        pass

    def minimise(self, stopping_dm=1e-2, max_steps=2000,
                 save_data_steps=10, save_m_steps=None, save_vtk_steps=None,
                 log_steps=1000):
        pass

    def relax(self):
        print('Not implemented for the SD minimiser')

    # -------------------------------------------------------------------------

    def compute_effective_field(self, t=0):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the corresponding Energy class
        """

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t=0, spin=self.spin)

    # -------------------------------------------------------------------------
    # SAVERS

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
        self.VTK.save_scalar(self._magnetisation / const.mu_B,
                             name='magnetisation')
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
            with zipfile.ZipFile('%s_m.zip' % self.name, 'a') as myzip:
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
