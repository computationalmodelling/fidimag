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
        self.nxyz = mesh.nxyz

        self.mu_s = mu_s
        self.mu_s_long = np.zeros(3 * mesh.nxyz)

        self.mu_s_long.shape = (-1, 3)
        for i in range(mesh.nxyz):
            self.mu_s_long[i, :] = mu_s[i]

        self.mu_s_long.shape = (-1,)

        self.field = np.zeros(3 * self.nxyz)
        self.field[:] = helper.init_vector(self.H0, self.mesh)

    def compute_field(self, t=0):
        return self.field

    # Todo: update it later
    def average_field(self):
        return self.field[0:3]

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
