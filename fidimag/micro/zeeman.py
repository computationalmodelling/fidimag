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

    # Todo: update it later
    def average_field(self):
        hx = self.field[0]
        hy = self.field[self.n]
        hz = self.field[2 * self.n]
        return np.array([hx, hy, hz])

    def compute_energy(self):

        sf = self.field * self.spin * self.Ms_long * mu_0

        energy = -np.sum(sf)

        return energy * self.mesh.cellsize


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
