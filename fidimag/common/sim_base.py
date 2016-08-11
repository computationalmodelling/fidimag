import os
import fidimag.common.helper as helper
from fidimag.common.fileio import DataSaver, DataReader

import numpy as np


class SimBase(object):
    """

    A class with common methods and definitions for both micromagnetic and
    atomistic simulations

    """
    def __init__(self, mesh, name):

        self.name = name
        self.mesh = mesh
        self.n = mesh.n
        self.n_nonzero = mesh.n
        self.unit_length = mesh.unit_length

        self.t = 0
        self.step = 0

        self._magnetisation = np.zeros(self.n, dtype=np.float)
        self.spin = np.ones(3 * self.n, dtype=np.float)
        self._pins = np.zeros(self.n, dtype=np.int32)
        self.field = np.zeros(3 * self.n, dtype=np.float)
        self._skx_number = np.zeros(self.n, dtype=np.float)
        self.interactions = []

        # This is for old C files codes using the xperiodic variables
        try:
            self.xperiodic, self.yperiodic, self.zperiodic = mesh.periodicity
        except ValueError:
            self.xperiodic, self.yperiodic = mesh.periodicity

        # To save the simulation data: ----------------------------------------

        self.data_saver = DataSaver(self, name + '.txt')

        self.data_saver.entities['E_total'] = {
            'unit': '<J>',
            'get': lambda sim: self.compute_energy(),
            'header': 'E_total'}

        self.data_saver.entities['m_error'] = {
            'unit': '<>',
            'get': lambda sim: self.compute_spin_error(),
            'header': 'm_error'}

        self.data_saver.update_entity_order()

        # ---------------------------------------------------------------------

    def set_m(self, m0=(1, 0, 0),
              normalise=True):
        """

        Set the magnetisation/spin three dimensional vector field profile.

        ARGUMENTS:

        m0      :: * To set every spin with the same direction,
                   set this value as a 3 elements tuple or list.

                   * For a spatially dependent vector field, you can specify a
                   function that returns a 3 element list depending on the
                   spatial coordinates. For example, a magnetisation field that
                   depends on the x position:

                        def m_profile(r):
                            for r[0] > 2:
                                return (0, 0, 1)
                            else:
                                return (0, 0, -1)

                   * You can also manually specify an array with (3 * n)
                   elements with the spins directions in the following order:

                        [mx_0 my_0 mz_0 mx_1 my_1 ... mx_n, my_n, mz_n]

                   where n is the number of mesh nodes and the order of the
                   magnetisation vectors follow the same order than the mesh
                   coordinates array.

                   * Alternatively, if you previously saved the magnetisation
                   field array to a numpy file, you can load it using
                   numpy.load(my_array)

        """

        self.spin[:] = helper.init_vector(m0, self.mesh, normalise)

        # TODO: carefully checking and requires to call set_mu first

        # Set the magnetisation/spin direction to (0, 0, 0) for sites
        # with no material, i.e. M_s = 0 or mu_s = 0
        # TODO: Check for atomistic and micromagnetic cases
        self.spin.shape = (-1, 3)
        for i in range(self.spin.shape[0]):
            if self._magnetisation[i] == 0:
                self.spin[i, :] = 0
        self.spin.shape = (-1,)

        # Set the initial state for the Sundials integrator using the
        # spins array
        # we probably need to change this line in the future.
        self.driver.set_initial_value(self.t)

    def get_pins(self):
        """
        Returns the array with pinned spins in the sample:
        sites with 0 are unpinned and sites with 1 are pinned. The order
        of the array follows the same order than the mesh coordinates
        """
        return self._pins

    def set_pins(self, pin):
        """

        An scalar field with values 1 or 0 to specify mesh/lattice sites
        with pinned or unpinned magnetic moments, respectively

        ARGUMENTS:

        pin     :: * You can specify a function that returns 1 or 0 depending
                   on the spatial coordinates. For example, to pin the spins
                   in a range in the x direction:

                        def pin_profile(r):
                            for r[0] > 2 and r[0] < 4:
                                return 1
                            else:
                                return 0

                   * You can also manually specify an array with n elements (1
                   or 0) with the pinned/unpinned values in the same order than
                   the mesh coordinates array.

                   * Alternatively, if you previously saved the pin
                   field array to a numpy file, you can load it using
                   numpy.load(my_array)
        """
        self._pins[:] = helper.init_scalar(pin, self.mesh)

        # Sites with no material, i.e. Mu_s or mu_s equal to zero,
        # will be pinned
        for i in range(len(self._magnetisation)):
            if self._magnetisation[i] == 0.0:
                self._pins[i] = 1

    pins = property(get_pins, set_pins)

    def get_alpha(self):
        """
        Returns the array with the spatially dependent Gilbert damping
        per mesh/lattice site
        """
        return self.driver._alpha

    def set_alpha(self, value):
        """

        Set the Gilbert damping of the system as a uniform or spatially
        dependent scalar field

        ARGUMENTS:

        value     :: * For a uniform damping across the whole sample, just
                       specify a float ranging from 0 to 1

                     * In addition, you can specify a function that returns
                       values ranging from 0 to 1, which depends on the spatial
                       coordinates. For example, a damping that increases
                       linearly in the x direction:

                        def alpha_profile(r):
                            for r[0] <= 10:
                                return r[0] / 10.
                            else:
                                return 0

                     * You can also manually specify an array with n values
                     ranging from 0 to 1 with the damping values, in the same
                     order than the mesh coordinates array.

                     * Alternatively, if you previously saved the damping
                     field array to a numpy file, you can load it using
                     numpy.load(my_array)
        """
        self.driver._alpha[:] = helper.init_scalar(value, self.mesh)
        #TODO: will encourge the users to set alpha through driver.alpha

    alpha = property(get_alpha, set_alpha)

    def add(self, interaction, save_field=False):
        """

        Add an interaction from one of the Energy subclasses. By default,
        the average energy of the added interaction is saved to the
        data file when relaxing the system

        OPTIONAL ARGUMENTS:

        save_field      :: Set True to save the average values of this
                           interaction field when relaxing the system

        """

        # magnetisation is Ms for the micromagnetic Sim class, and it is
        # mu_s for the atomistic Sim class
        interaction.setup(self.mesh, self.spin,
                          self._magnetisation
                          )

        # TODO: FIX  --> ??
        # When adding an interaction that was previously added, using
        # the same name, append a '_2' to the new interaction name (?)
        for i in self.interactions:
            if i.name == interaction.name:
                interaction.name = i.name + '_2'

        self.interactions.append(interaction)

        # Specify a name for the energy of the interaction, which will
        # appear in a file with saved values
        # When saving the energy values, we call the compute_energy() method
        # from the (micromagnetic/atomistic) Energy class (overhead?)
        energy_name = 'E_{0}'.format(interaction.name)
        self.data_saver.entities[energy_name] = {
            'unit': '<J>',
            'get': lambda sim: sim.get_interaction(interaction.name).compute_energy(),
            'header': energy_name}

        # Save the average values of the interaction vector field components
        if save_field:
            fn = '{0}'.format(interaction.name)
            self.data_saver.entities[fn] = {
                'unit': '<>',
                'get': lambda sim: sim.get_interaction(interaction.name).average_field(),
                'header': ('%s_x' % fn, '%s_y' % fn, '%s_z' % fn)}

        self.data_saver.update_entity_order()

    def get_interaction(self, name):
        """
        Returns an instance of a magnetic interaction previously added
        to the simulation, using the corresponding interaction name as
        a string
        """
        for interaction in self.interactions:
            if interaction.name == name:
                return interaction
        else:
            raise ValueError("Failed to find the interaction with name '{0}', "
                             "available interactions: {1}.".format(
                                 name, [x.name for x in self.interactions]))

    def skyrmion_number(self):
        pass

    def spin_at(self, i, j, k):
        """
        Returns the x,y,z components of a spin in the [i, j, k]
        position of the mesh, where i,j,k are integer indexes. The index
        ordering is specified in the mesh class.
        """

        i1 = 3 * self.mesh.index(i, j, k)

        # print self.spin.shape,nxy,nx,i1,i2,i3
        return np.array([self.spin[i1],
                         self.spin[i1 + 1],
                         self.spin[i1 + 2]])

    def add_monitor_at(self, i, j, k, name='p'):
        """
        Save site spin with index (i,j,k) to txt file.
        """

        self.data_saver.entities[name] = {
            'unit': '<>',
            'get': lambda sim: sim.spin_at(i, j, k),
            'header': (name + '_x', name + '_y', name + '_z')}

        self.data_saver.update_entity_order()

    def spin_length(self):
        """
        Returns an array with the length of every spin in the mesh. The
        order is given by the mesh.coordinates order
        """
        self.spin.shape = (-1, 3)
        length = np.sqrt(np.sum(self.spin ** 2, axis=1))
        self.spin.shape = (-1,)
        return length

    def compute_spin_error(self):
        length = self.spin_length() - 1.0
        length[self._pins > 0] = 0
        return np.max(abs(length))

    def compute_average(self):
        """
        Compute the average values of the 3 components of the magnetisation
        vector field
        """
        self.spin.shape = (-1, 3)
        average = np.sum(self.spin, axis=0) / self.n_nonzero
        self.spin.shape = (3 * self.n)
        return average

    def compute_energy(self):
        """
        Compute the total energy of the magnetic system
        """
        energy = 0

        for obj in self.interactions:
            energy += obj.compute_energy()

        return energy

    # -------------------------------------------------------------------------
    # Save functions ----------------------------------------------------------
    # -------------------------------------------------------------------------

    def save_vtk(self):
        pass

    def save_m(self):
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

    def run_until(self, t):
        """
        Evolve the system with a micromagnetic driver (LLG, LLG_STT, etc.)
        until a specific time `t`, using the specified integrator.
        The integrator was specified with the right hand side of the
        driver equation

        """

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
            return

        self.driver.next_step(t)

        self.t = self.driver.t
        self.step = self.driver.step

        # Update field before saving data
        self.compute_effective_field(t)
        self.data_saver.save()

    def relax(self, dt=1e-11, stopping_dmdt=0.01, max_steps=1000,
              save_m_steps=100, save_vtk_steps=100
              ):
        """

        Evolve the system until meeting the dmdt < stopping_dmdt criteria. We
        can specify to save VTK and NPY files with the magnetisation vector
        field, every certain number of the integrator steps

        TODO: Check what dt is exactly doing

        """

        for i in range(0, max_steps + 1):

            cvode_dt = self.driver.cvode_dt()

            increment_dt = dt
            if dt < cvode_dt:
                increment_dt = cvode_dt
                self.driver.next_step()
            else:
                self.driver.next_step(self.t+dt)
            
            self.t = self.driver.t
            self.step = self.driver.step

            if save_vtk_steps is not None:
                if i % save_vtk_steps == 0:
                    self.save_vtk()
            if save_m_steps is not None:
                if i % save_m_steps == 0:
                    self.save_m()

            dmdt = self.driver.compute_dmdt(increment_dt)

            print(('step=%d ' +
                   'time=%0.3g ' +
                   'max_dmdt=%0.3g ' +
                   'ode_step=%0.3g') % (self.step,
                                        self.t,
                                        dmdt / self.driver._dmdt_factor,
                                        cvode_dt)
                  )


            if dmdt < stopping_dmdt * self.driver._dmdt_factor:
                break

        if save_m_steps is not None:
            self.save_m()

        if save_vtk_steps is not None:
            self.save_vtk()
