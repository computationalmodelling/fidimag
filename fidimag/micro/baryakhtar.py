import numpy as np
import fidimag.extensions.baryakhtar_clib as clib
from .micro_driver import MicroDriver
from .relax import Relaxation
from .relax import Laplace


class LLBarFull(MicroDriver):

    def __init__(self, mesh, spin, Ms, Ms_inv, field, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False,
                 chi=1e-3
                 ):

        # Inherit from the driver class
        super(LLBarFull, self).__init__(mesh, spin, Ms, Ms_inv, field,
                                        pins, interactions, name,
                                        data_saver,
                                        integrator=integrator,
                                        use_jac=use_jac
                                        )

        self._chi = 1e-3
        add(Relaxation(self.chi, name='chi_relax'),
            self.interactions, self.data_saver,
            self.mesh, self.spin, self._Ms)

        self.lap = Laplace(mesh)
        # OLD: self.add(Relaxation(chi))

        self.beta = 0

    def get_chi(self):
        """
        """
        return self._chi

    def set_chi(self, chi):
        """
        """
        for i in self.interactions:
            if i.name == 'chi_relax':
                i.chi = chi

    chi = property(get_chi, set_chi)

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)
        delta_h = self.lap.compute_laplace_field(self.field, self._Ms)

        clib.compute_llg_rhs_baryakhtar(ydot,
                                        self.spin,
                                        self.field,
                                        delta_h,
                                        self.alpha,
                                        self.beta,
                                        self._pins,
                                        self.gamma,
                                        self.n,
                                        self.do_precession)

        #ydot[:] = self.dm_dt[:]

        return 0


class LLBar(MicroDriver):

    def __init__(self, mesh, spin, Ms, field, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False,
                 ):

        # Inherit from the driver class
        super(LLBar, self).__init__(mesh, spin, Ms, field,
                                    pins, interactions, name,
                                    data_saver,
                                    integrator='sundials',
                                    use_jac=False
                                    )

        self.lap = Laplace(mesh)

        self.field_perp = np.zeros(3 * self.n, dtype=np.float64)

        self.beta = 0

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)
        clib.compute_perp_field(
            self.spin, self.field, self.field_perp, self.n)
        delta_h = self.lap.compute_laplace_field(self.field_perp, self._Ms)

        clib.compute_llg_rhs_baryakhtar_reduced(ydot,
                                                self.spin,
                                                self.field,
                                                delta_h,
                                                self.alpha,
                                                self.beta,
                                                self._pins,
                                                self.gamma,
                                                self.n,
                                                self.do_precession,
                                                self.default_c)

        #ydot[:] = self.dm_dt[:]

        return 0


def add(interaction, interactions_list, data_saver,
        mesh, spin, magnetisation, save_field=False):

    """

    Add an interaction to the interaction list.
    This function is based on the Sim class *add* method
    (it is likely that this function will be moved to the common
     helpers in the future)

    OPTIONAL ARGUMENTS:

    save_field      :: Set True to save the average values of this
                       interaction field when relaxing the system

    """

    # magnetisation is Ms for the micromagnetic Sim class, and it is
    # mu_s for the atomistic Sim class
    interaction.setup(mesh, spin, magnetisation)

    # TODO: FIX  --> ??
    # When adding an interaction that was previously added, using
    # the same name, append a '_2' to the new interaction name (?)
    for i in interactions_list:
        if i.name == interaction.name:
            interaction.name = i.name + '_2'

    interactions_list.append(interaction)

    # Specify a name for the energy of the interaction, which will
    # appear in a file with saved values
    # When saving the energy values, we call the compute_energy() method
    # from the (micromagnetic/atomistic) Energy class (overhead?)
    energy_name = 'E_{0}'.format(interaction.name)
    data_saver.entities[energy_name] = {
        'unit': '<J>',
        'get': lambda sim: sim.get_interaction(interaction.name).compute_energy(),
        'header': energy_name}

    # Save the average values of the interaction vector field components
    if save_field:
        fn = '{0}'.format(interaction.name)
        data_saver.entities[fn] = {
            'unit': '<>',
            'get': lambda sim: sim.get_interaction(interaction.name).average_field(),
            'header': ('%s_x' % fn, '%s_y' % fn, '%s_z' % fn)}

    data_saver.update_entity_order()
