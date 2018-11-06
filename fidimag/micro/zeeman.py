import numpy as np
from fidimag.common.constant import mu_0
import fidimag.common.helper as helper
import inspect


class Zeeman(object):

    """
    A time independent external magnetic field that can be space dependent.
    The field energy in the micromagnetic theory reads:
                                 _
                               /   ->       ->
           E   =  - \mu_0     /    M  \cdot H   dV
                           _ /

    with H as the bias field in A / m, \mu_0 the vacuum permeability and M the
    magnetisation vector. Using finite differences, this quantity is computed
    through the summation

                              __  ->        ->
         E =  - \mu_0 * dV   \    M_i \cdot H_i
                             /__
                              i

    where M_i is the magnetisation at the i-th position of the mesh
    discretisation and dV = dx * dy * dz is the volume of a mesh unit cell.

    If the field is homogeneous, it can be specified in a simulation object
    *Sim* as

            Sim.add(Zeeman((H_x, H_y, H_z)))

    Otherwise, it can be specified as any Fidimag field, passing a function or
    an array. For example, a space dependent field function that changes
    linearly in the x-direction, and only has a x-component, can be defined as:

        def my_Zeeman_field(pos):
            H = 0.01 / (4 * np.pi * 1e-7)  # A / m
            return (H * pos[0], 0, 0)

        # Add field to Simulation object
        Sim.add(Zeeman(my_Zeeman_field))

    For a hysteresis loop, the field can be updated using the *update_field*
    function. For an already defined simulation object *Sim*, it is updated as

        Sim.get_interaction('Zeeman').update_field(new_field)

    where new_field is a 3-tuple, and array or a function, as shown before.
    It is recommended to reset the integrator with *Sim.reset_integrator()*
    to start a new relaxation after updating the field.

    """

    def __init__(self, H0, name='Zeeman'):
        # Raise an exception if H0 is not indexable or callable. This is
        # because H0 represents a vector.
        if hasattr(H0, "__getitem__") is False and\
           hasattr(H0, "__call__") is False:
            raise ValueError("H0 \"{}\" does not represent a vector"
                             .format(H0))
        self.H0 = H0
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, Ms, Ms_inv):
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n

        self.Ms = Ms

        # TODO: Check if it is necessary to define a 3D matrix for
        # the Ms vectors. Maybe there is a way that uses less memory
        # (see the calculation in the *compute_energy* function)
        self.field = np.zeros(3 * self.n)
        self.field[:] = helper.init_vector(self.H0, self.mesh, 3)
        # print self.field

    def update_field(self, H0):
        self.H0 = H0
        self.field[:] = helper.init_vector(self.H0, self.mesh, 3)

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

        sf = self.field * self.spin * mu_0

        energy_density = -np.sum(sf.reshape(-1, 3), axis=1) * self.Ms

        return np.sum(energy_density) * (self.mesh.dx *
                                         self.mesh.dy *
                                         self.mesh.dz *
                                         self.mesh.unit_length ** 3.)


class TimeZeeman(Zeeman):

    """
    The time dependent external field, also can vary with space
    """

    def __init__(self, H0, time_fun, extra_args=[], name='TimeZeeman'):
        self.H0 = H0
        self.time_fun = time_fun
        self.name = name
        self.jac = True
        self.extra_args = extra_args

    def setup(self, mesh, spin, Ms, Ms_inv):
        super(TimeZeeman, self).setup(mesh, spin, Ms, Ms_inv)
        self.H_init = self.field.copy()

    def compute_field(self, t=0, spin=None):
        self.field[:] = self.H_init[:] * self.time_fun(t, *self.extra_args)
        return self.field
