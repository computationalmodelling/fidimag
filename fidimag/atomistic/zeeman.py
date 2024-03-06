import numpy as np
import fidimag.common.helper as helper


class Zeeman(object):

    """

    A time independent external magnetic field that can be space dependent.
    The field energy is computed as:

                  __   ->       ->
         E =  -  \    \mu_i  *  B_i
                 /__
                  i

    where mu_i = g \mu_B S_i is the magnetic moment vector at the i-th lattice
    site, g is the Lande factor, \mu_B the Bohr magneton, S_i the average total
    spin vector at the i-th site and B_i the bias field vector at the i-th
    site, given in Tesla units.

    OPTIONAL ARGUMENTS: -------------------------------------------------------
        name            :: Interaction name

    USAGE: --------------------------------------------------------------------

    If the field is homogeneous, it can be specified in a simulation object
    *Sim* as

            Sim.add(Zeeman((B_x, B_y, B_z)))

    Otherwise, it can be specified as any Fidimag vector field, passing a
    function or an array. For example, a space dependent field function that
    changes linearly in the x-direction, and only has a x-component, can be
    defined as:

        def my_Zeeman_field(pos):
            B = 0.01  # T
            return (B * pos[0], 0, 0)

        # Add field to Simulation object
        Sim.add(Zeeman(my_Zeeman_field))


    For a hysteresis loop, the field can be updated using the *update_field*
    function. For an already defined simulation object *Sim*, it is updated as

        Sim.get_interaction('Zeeman').update_field(new_field)

    where new_field is a 3-tuple, and array or a function, as shown before.
    It is recommended to reset the integrator with *Sim.reset_integrator()*
    to start a new relaxation after updating the field.

    """

    def __init__(self, B0, name='Zeeman'):
        self.B0 = B0
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n

        self.mu_s = mu_s

        self.field = np.zeros(3 * self.n)
        self.field[:] = helper.init_vector(self.B0, self.mesh)
        self.energy = np.zeros(mesh.n)

    def update_field(self, B0):
        self.B0 = B0
        self.field[:] = helper.init_vector(self.B0, self.mesh)

    def compute_field(self, t=0, spin=None):
        return self.field

    def average_field(self):
        # Remember that fields are: [fx0, fy0, fz0, fx1, fy1, fz1, fx2, ...]
        # So we jump in steps of 3 starting from the 0, 1 and 2nd elements
        hx = self.field[::3]
        hy = self.field[1::3]
        hz = self.field[2::3]
        return np.array([np.average(hx), np.average(hy), np.average(hz)])

    def compute_energy(self):

        sf = self.field * self.spin
        energy_density = -np.sum(sf.reshape(-1, 3), axis=1) * self.mu_s
        self.energy[:] = energy_density
        self.total_energy = np.sum(self.energy)

        return self.total_energy


class TimeZeeman(Zeeman):

    """
    The time dependent external field, also can vary with space
    """

    def __init__(self, B0, time_fun, name='TimeZeeman'):
        self.B0 = B0
        self.time_fun = time_fun
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(TimeZeeman, self).setup(mesh, spin, mu_s, mu_s_inv)
        self.H_init = self.field.copy()

    def compute_field(self, t=0):
        self.field[:] = self.H_init[:] * self.time_fun(t)
        return self.field
