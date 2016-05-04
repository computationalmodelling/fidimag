import numpy as np
import fidimag.common.helper as helper


class Zeeman(object):

    """
    The time independent external field, can vary with space
    """

    def __init__(self, H0, name='Zeeman'):
        self.H0 = H0
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, mu_s):
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n

        self.mu_s = mu_s
        self.mu_s_long = np.zeros(3 * mesh.n)

        self.mu_s_long.shape = (-1, 3)
        for i in range(mesh.n):
            self.mu_s_long[i, :] = mu_s[i]

        self.mu_s_long.shape = (-1,)

        self.field = np.zeros(3 * self.n)
        self.field[:] = helper.init_vector(self.H0, self.mesh)

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

        sf = self.field * self.spin * self.mu_s_long

        energy = -np.sum(sf)

        return energy


class TimeZeeman(Zeeman):

    """
    The time dependent external field, also can vary with space
    """

    def __init__(self, H0, time_fun, name='TimeZeeman'):
        self.H0 = H0
        self.time_fun = time_fun
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, mu_s):
        super(TimeZeeman, self).setup(mesh, spin, mu_s)
        self.H_init = self.field.copy()

    def compute_field(self, t=0):
        self.field[:] = self.H_init[:] * self.time_fun(t)
        return self.field
