import numpy as np

from fidimag.common.constant import mu_0
import fidimag.common.helper as helper


class Zeeman(object):

    """
    The time independent external field, can vary with space
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

    def setup(self, mesh, spin, Ms):
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n

        self.Ms = Ms
        self.Ms_long = np.zeros(3 * mesh.n)

        self.Ms_long.shape = (3, -1)
        for i in range(mesh.n):
            self.Ms_long[:, i] = Ms[i]

        self.Ms_long.shape = (-1,)
        self.field = np.zeros(3 * self.n)
        self.field[:] = helper.init_vector(self.H0, self.mesh)
        # print self.field

    def update_field(self, H0):
        self.H0 = H0
        self.field[:] = helper.init_vector(self.H0, self.mesh)

    def compute_field(self, t=0):
        return self.field

    def average_field(self):
        # Remember that fields are: [fx0, fy0, fz0, fx1, fy1, fz1, fx2, ...]
        # So we jump in steps of 3 starting from the 0, 1 and 2nd elements
        hx = self.field[::3]
        hy = self.field[1::3]
        hz = self.field[2::3]
        return np.array([np.average(hx), np.average(hy), np.average(hz)])

    def compute_energy(self):

        sf = self.field * self.spin * self.Ms_long * mu_0

        energy = -np.sum(sf)

        return energy * (self.mesh.dx *
                         self.mesh.dy *
                         self.mesh.dz *
                         self.mesh.unit_length ** 3.)


class TimeZeeman(Zeeman):

    """
    The time dependent external field, also can vary with space
    """

    def __init__(self, H0, time_fun, name='TimeZeeman'):
        self.H0 = H0
        self.time_fun = time_fun
        self.name = name

    def setup(self, mesh, spin, Ms):
        super(TimeZeeman, self).setup(mesh, spin, Ms)
        self.H_init = self.field.copy()

    def compute_field(self, t=0):
        self.field[:] = self.H_init[:] * self.time_fun(t)
        return self.field
