import numpy as np


class Energy:

    """
    An abstract class to implement the basic functions such as setup in micromagnetics.
    """

    def setup(self, mesh, spin, Ms, Ms_inv):
        self.mesh = mesh
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.dxyz = self.dx * self.dy * self.dz
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

    # WARNING:this should be time-dependent in compute_field
    def compute_energy(self):
        """
        Energy of the whole sample, and of every cell

        Fills `self.energy` with the energy of each cell, in joules, and sets
        `self.total_energy` to their sum. The two are kept in that relation by
        every interaction class, micromagnetic and atomistic, so that a per
        cell energy can be summed or differenced without knowing which
        interaction produced it.

        The `compute_field` routines write an energy *density* into
        `self.energy`, since that is what the finite difference expressions
        give, so the cell volume is applied here.

        Returns
        -------
        float
            `self.total_energy`
        """
        # since we are not always calling this function, so it's okay to call
        # compute_field again
        self.compute_field()

        self.energy *= self.dxyz
        self.total_energy = np.sum(self.energy)

        return self.total_energy
