import numpy as np


class Energy(object):

    """
    An abstract class to implement the basic functions such as setup in micromagnetics.
    """

    def setup(self, mesh, spin, Ms, Ms_inv):
        self.mesh = mesh
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = spin
        self.n = mesh.n

        self.field = np.zeros(3 * mesh.n)
        self.energy = np.zeros(mesh.n)
        self.total_energy = 0
        self.Ms = Ms
        self.Ms_inv = Ms_inv

        # For old code compatibility
        self.xperiodic, self.yperiodic, self.zperiodic = mesh.periodicity

        self.neighbours = mesh.neighbours

    def compute_field(self, t=0):

        return 0

    def compute_energy(self):

        # since we are not always calling this function, so it's okay to call
        # compute_field again
        self.compute_field()

        self.total_energy = np.sum(self.energy) * (self.mesh.dx *
                                                   self.mesh.dy *
                                                   self.mesh.dz *
                                                   self.mesh.unit_length ** 3.)

        return self.total_energy
