from . import llg
from . import llg_stt
from . import llg_stt_cpp
from . import baryakhtar

from . import micro_driver

import numpy as np

KNOWN_DRIVERS = {'llg': llg.LLG,
                 'llg_stt': llg_stt.LLG_STT,
                 'llg_stt_cpp': llg_stt_cpp.LLG_STT_CPP,
                 'llbar': baryakhtar.LLBar,
                 'llbar_full': baryakhtar.LLBarFull}


class Sim(object):
    """
    def Sim(*args, **kwargs)
    
    This function returns a simulation object.

    By default, it will return a simulation object for the LLG driver;
    other available drivers are:
        llg_stt - LLG w. spin transfer torque
        llg_stt_cpp - LLG w. spin transfer torque 
        llbar - Landau-Lifshitz-Baryakhtar equation
        llbar_full
    
    To use the STT driver, for example, pass as an argument driver="llg_stt".

    """
    def __init__(self, mesh,
                 driver='llg', 
                 *args, **kwargs
                 ):
        self._micromagnetic = True

        self.name = name
        self.mesh = mesh
        self.n = mesh.n
        self.n_nonzero = mesh.n
        self.unit_length = mesh.unit_length


        self._Ms = np.zeros(self.n, dtype=np.float)
        self.spin = np.ones(3 * self.n, dtype=np.float)
        self._pins = np.zeros(self.n, dtype=np.int32)
        self.field = np.zeros(3 * self.n, dtype=np.float)

        self._skx_number = np.zeros(self.n, dtype=np.float)

        # This is for old C files codes using the xperiodic variables
        self.xperiodic, self.yperiodic, self.zperiodic = mesh.periodicity

        # 
        self.driver = KNOWN_DRIVERS[driver](*args, **kwargs)


    def set_m(self, m0=(1, 0, 0), normalise=True):

        self.spin[:] = helper.init_vector(m0, self.mesh, normalise)

        # TODO: carefully checking and requires to call set_mu first
        self.spin.shape = (-1, 3)
        for i in range(self.spin.shape[0]):
            if self._Ms[i] == 0:
                self.spin[i, :] = 0
        self.spin.shape = (-1,)

        self.integrator.set_initial_value(self.spin, self.t)

    def get_pins(self):
        return self._pins

    def set_pins(self, pin):
        self._pins[:] = helper.init_scalar(pin, self.mesh)

        for i in range(len(self._Ms)):
            if self._Ms[i] == 0.0:
                self._pins[i] = 1

    pins = property(get_pins, set_pins)

    def get_alpha(self):
        return self._alpha

    def set_alpha(self, value):
        self._alpha[:] = helper.init_scalar(value, self.mesh)

    alpha = property(get_alpha, set_alpha)

    def get_Ms(self):
        return self._Ms

    def set_Ms(self, value):
        self._Ms[:] = helper.init_scalar(value, self.mesh)
        nonzero = 0
        for i in range(self.n):
            if self._Ms[i] > 0.0:
                self._Ms_inv = 1.0 / self._Ms[i]
                nonzero += 1

        self.n_nonzero = nonzero

        for i in range(len(self._Ms)):
            if self._Ms[i] == 0.0:
                self._pins[i] = 1

                # Set the neighbour index to -1 for sites with Ms = 0
                self.mesh.neighbours[self.mesh.neighbours == i] = -1

        self.Ms_const = np.max(self._Ms)

    Ms = property(get_Ms, set_Ms)

    def add(self, interaction, save_field=False):
        interaction.setup(self.mesh, self.spin, Ms=self._Ms)

        # TODO: FIX
        for i in self.interactions:
            if i.name == interaction.name:
                interaction.name = i.name + '_2'

        self.interactions.append(interaction)

        energy_name = 'E_{0}'.format(interaction.name)
        self.saver.entities[energy_name] = {
            'unit': '<J>',
            'get': lambda sim: sim.get_interaction(interaction.name).compute_energy(),
            'header': energy_name}

        if save_field:
            fn = '{0}'.format(interaction.name)
            self.saver.entities[fn] = {
                'unit': '<>',
                'get': lambda sim: sim.get_interaction(interaction.name).average_field(),
                'header': ('%s_x' % fn, '%s_y' % fn, '%s_z' % fn)}

        self.saver.update_entity_order()

    def get_interaction(self, name):
        for interaction in self.interactions:
            if interaction.name == name:
                return interaction
        else:
            raise ValueError("Failed to find the interaction with name '{0}', "
                             "available interactions: {1}.".format(
                                 name, [x.name for x in self.interactions]))

    def compute_average(self):
        self.spin.shape = (-1, 3)
        average = np.sum(self.spin, axis=0) / self.n_nonzero
        self.spin.shape = (3 * self.n)
        return average

    def compute_energy(self):

        energy = 0

        for obj in self.interactions:
            energy += obj.compute_energy()

        return energy

    def skyrmion_number(self):
        nx = self.mesh.nx
        ny = self.mesh.ny
        nz = self.mesh.nz
        number = micro_clib.compute_skyrmion_number(
            self.spin, self._skx_number, nx, ny, nz, self.mesh.neighbours)
        return number

    skyrmion_number_slice = fidimag.common.skyrmion_number.skyrmion_number_slice
    skyrmion_number_lee = fidimag.common.skyrmion_number.skyrmion_number_lee

    def spin_at(self, i, j, k):

        i1 = 3 * self.mesh.index(i, j, k)

        # print self.spin.shape,nxy,nx,i1,i2,i3
        return np.array([self.spin[i1],
                         self.spin[i1 + 1],
                         self.spin[i1 + 2]])

    def add_monitor_at(self, i, j, k, name='p'):
        """
        Save site spin with index (i,j,k) to txt file.
        """

        self.saver.entities[name] = {
            'unit': '<>',
            'get': lambda sim: sim.spin_at(i, j, k),
            'header': (name + '_x', name + '_y', name + '_z')}

        self.saver.update_entity_order()

    def save_vtk(self):
        self.vtk.save_vtk(self.spin.reshape(-1, 3), self.Ms, step=self.step)

    def save_m(self):
        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)

    def save_skx(self):
        if not os.path.exists('%s_skx_npys' % self.name):
            os.makedirs('%s_skx_npys' % self.name)
        name = '%s_skx_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self._skx_number)

    def stat(self):
        return self.integrator.stat()

    def spin_length(self):
        self.spin.shape = (-1, 3)
        length = np.sqrt(np.sum(self.spin**2, axis=1))
        self.spin.shape = (-1,)
        return length

    def compute_spin_error(self):
        length = self.spin_length() - 1.0
        length[self._pins > 0] = 0
        return np.max(abs(length))

    def compute_dmdt(self, dt):
        m0 = self.spin_last
        m1 = self.spin
        dm = (m1 - m0).reshape((3, -1))
        max_dm = np.max(np.sqrt(np.sum(dm**2, axis=0)))
        max_dmdt = max_dm / dt
        return max_dmdt

def Sim(*args, **kwargs):

    driver = 'llg'

    if 'driver' in kwargs:
        driver = kwargs.pop('driver')

    if driver not in KNOWN_DRIVERS:
        raise NotImplementedError("""Driver '{}' is not implemented.
                                  Valid choices: one of '{}'.""".format(driver, KNOWN_DRIVERS.keys()))

    return KNOWN_DRIVERS[driver](*args, **kwargs)
