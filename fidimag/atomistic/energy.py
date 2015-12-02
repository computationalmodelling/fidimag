import numpy as np

#from constant import mu_0


class Energy(object):

    """

    An abstract class to implement the basic functions such as setup for atomic
    dynamics.

    """

    def setup(self, mesh, spin, mu_s):
        self.mesh = mesh
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.n = mesh.n
        self.mesh_type = mesh.mesh_type

        self.total_energy = 0
        self.mu_s = mu_s

        self.field = np.zeros(3 * self.n, dtype=np.float)
        self.energy = np.zeros(3 * mesh.n, dtype=np.float)
        self.mu_s_inv = np.zeros(3 * self.n, dtype=np.float)

        self.mu_s_inv.shape = (-1, 3)
        for i in range(mesh.n):
            if self.mu_s[i] == 0.0:
                self.mu_s_inv[i, :] = 0
            else:
                self.mu_s_inv[i, :] = 1.0 / self.mu_s[i]

        self.mu_s_inv.shape = (-1,)

        self.xperiodic, self.yperiodic = (mesh.periodicity[0],
                                          mesh.periodicity[1])
        self.neighbours = mesh.neighbours

        try:
            self.coordinates = self.mesh.coordinates
        except:
            self.coordinates = np.array(self.mesh.pos)

    def compute_field(self, t=0):

        return 0

    def compute_energy(self):

        # since we are not always calling this function, so it's okay to call
        # compute_field again
        self.compute_field()

        self.total_energy = np.sum(self.energy)

        return self.total_energy
