import fidimag.extensions.dipolar as clib
import numpy as np


class Demag(object):

    def __init__(self, name='demag'):
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s):
        self.mesh = mesh
        self.dx = mesh.dx
        self.dy = mesh.dy
        self.dz = mesh.dz
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.n = mesh.n
        self.field = np.zeros(3 * self.n, dtype=np.float)
        unit_length = mesh.unit_length
        self.mu_s_scale = np.zeros(mesh.n, dtype=np.float)

        # note that the 1e-7 comes from \frac{\mu_0}{4\pi}
        self.scale = 1e-7 / unit_length**3

        # could be wrong, needs carefully tests!!!
        self.mu_s_scale = mu_s * self.scale

        self.demag = clib.FFTDemag(self.dx, self.dy, self.dz,
                                   self.nx, self.ny, self.nz,
                                   False)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        self.demag.compute_field(m, self.mu_s_scale, self.field)
        return self.field

    def compute_exact(self):
        field = np.zeros(3 * self.n)
        self.demag.compute_exact(self.spin, self.mu_s_scale, field)
        return field

    def compute_energy(self):

        energy = self.demag.compute_energy(
            self.spin, self.mu_s_scale, self.field)

        return energy / self.scale
