from __future__ import division
from __future__ import print_function
import os
import numpy as np
import fidimag.extensions.a_clib as clib
from fidimag.common.fileio import DataSaver
import fidimag.common.helper as helper
import fidimag.common.constant as const
from fidimag.common.save_vtk import SaveVTK

class MonteCarlo(object):

    def __init__(self, mesh, name='unnamed'):
        self.mesh = mesh
        self.name = name
        self.n = mesh.n
        self.n_nonzero = self.n

        self.ngbs = mesh.neighbours

        self._mu_s = np.zeros(self.n, dtype=np.float)
        self.spin = np.ones(3 * self.n, dtype=np.float)
        self.spin_last = np.ones(3 * self.n, dtype=np.float)

        self.random_spin = np.zeros(3 * self.n, dtype=np.float)
        self._H = np.zeros(3 * self.n, dtype=np.float)

        self._skx_number = np.zeros(self.n, dtype=np.float)
        self.interactions = []

        self.create_tablewriter()
        self.vtk = SaveVTK(self.mesh, name=name)

        self.hexagonal_mesh = False
        if mesh.mesh_type == 'hexagonal':
            self.hexagonal_mesh =  True
            #FIX ME !!!!
            self.nngbs = np.copy(mesh.neighbours)
        else:
            self.nngbs = mesh.next_neighbours

        self.step = 0
        self.skx_num = 0
        self.mc = clib.monte_carlo()
        self.set_options()

    def set_options(self, J=50.0*const.k_B, J1=0, D=0, D1=0, Kc=0, H=None, seed=100, T=10.0, S=1):
        """
        J, D and Kc in units of Joule
        H in units of Tesla.
        S is the spin length
        """
        self.mc.set_seed(seed)
        self.J = J/const.k_B
        self.J1 = J1/const.k_B
        self.D = D/const.k_B
        self.D1 = D1/const.k_B
        self.T = T
        self.Kc = Kc/const.k_B
        self.mu_s =  1.0
        if H is not None:
            self._H[:] = helper.init_vector(H, self.mesh)
            self._H[:] = self._H[:]*const.mu_s_1*S/const.k_B

    def create_tablewriter(self):

        entities = {
            'step': {'unit': '<>',
                     'get': lambda sim: sim.step,
                     'header': 'step'},
            'm': {'unit': '<>',
                  'get': lambda sim: sim.compute_average(),
                  'header': ('m_x', 'm_y', 'm_z')},
            'skx_num':{'unit': '<>',
                      'get': lambda sim: sim.skyrmion_number(),
                      'header': 'skx_num'}
        }

        self.saver = DataSaver(self, self.name + '.txt', entities=entities)

        self.saver.update_entity_order()

    def set_m(self, m0=(1, 0, 0), normalise=True):
        self.spin[:] = helper.init_vector(m0, self.mesh, norm=normalise)

        # TODO: carefully checking and requires to call set_mu first

        self.spin.shape = (-1, 3)
        for i in range(self.spin.shape[0]):
            if self._mu_s[i] == 0:
                self.spin[i, :] = 0
        self.spin.shape = (-1,)

    def get_mu_s(self):
        return self._mu_s

    def set_mu_s(self, value):
        self._mu_s[:] = helper.init_scalar(value, self.mesh)
        nonzero = 0
        for i in range(self.n):
            if self._mu_s[i] > 0.0:
                nonzero += 1

        self.n_nonzero = nonzero

    mu_s = property(get_mu_s, set_mu_s)

    def compute_average(self):
        self.spin.shape = (-1, 3)
        average = np.sum(self.spin, axis=0) / self.n_nonzero
        self.spin.shape = (-1,)
        return average

    def skyrmion_number(self):
        nx = self.mesh.nx
        ny = self.mesh.ny
        nz = self.mesh.nz
        number = clib.compute_skyrmion_number(
            self.spin, self._skx_number, nx, ny, nz, self.mesh.neighbours, self.mesh.n_ngbs)
        self.skx_num = number
        return number


    def save_vtk(self):
        """
        Save a VTK file with the magnetisation vector field and magnetic
        moments as cell data. Magnetic moments are saved in units of
        Bohr magnetons

        NOTE: It is recommended to use a *cell to point data* filter in
        Paraview or Mayavi to plot the vector field
        """
        self.vtk.save_vtk(self.spin.reshape(-1, 3),
                          self._mu_s,
                          step=self.step)

    def save_m(self):
        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)


    def run(self, steps=1000, save_m_steps=100, save_vtk_steps=100, save_data_steps=1):

        if save_m_steps is not None:
            self.save_m()

        if save_vtk_steps is not None:
            self.save_vtk()

        for step in range(1, steps + 1):
            self.step += 1
            self.mc.run_step(self.spin, self.random_spin,
                             self.ngbs, self.nngbs, self.mesh.n_ngbs,
                             self.J, self.J1, self.D, self.D1, self._H, self.Kc,
                             self.n, self.T, self.hexagonal_mesh)
            if save_data_steps is not None:
                if step % save_data_steps == 0:
                    self.saver.save()
                    print("step=%d, skyrmion number=%0.9g"%(self.step, self.skx_num))

            if save_vtk_steps is not None:
                if step % save_vtk_steps == 0:
                    self.save_vtk()
            if save_m_steps is not None:
                if step % save_m_steps == 0:
                    self.save_m()




if __name__ == '__main__':
    pass
