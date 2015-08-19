import numpy as np


class Energy(object):

    """
    An abstract class to implement the basic functions such as setup in micromagnetics.
    """

    def setup(self, mesh, spin, Ms):
        self.mesh = mesh
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.nxyz = mesh.nxyz

        self.field = np.zeros(3 * mesh.nxyz)
        self.energy = np.zeros(mesh.nxyz)
        self.total_energy = 0
        self.pbc = mesh.pbc
        self.Ms = Ms
        self.Ms_inv = np.zeros(mesh.nxyz)

        for i in range(mesh.nxyz):
            if self.Ms[i] == 0.0:
                self.Ms_inv[i] = 0
            else:
                self.Ms_inv[i] = 1.0 / self.Ms[i]

        self.xperiodic = mesh.xperiodic
        self.yperiodic = mesh.yperiodic
        self.zperiodic = mesh.zperiodic

    def compute_field(self, t=0):

        return 0

    def compute_energy(self):

        # since we are not always calling this function, so it's okay to call
        # compute_field again
        self.compute_field()

        self.total_energy = np.sum(self.energy) * self.mesh.cellsize

        return self.total_energy
