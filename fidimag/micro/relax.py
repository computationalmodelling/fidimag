import fidimag.extensions.baryakhtar_clib as clib
import numpy as np


class Relaxation:

    """
    compute the relaxation field related to exchange field.
    """

    def __init__(self, chi, name='relax'):
        self.chi = chi
        self.name = name

    def setup(self, mesh, m, Ms):
        self.mesh = mesh
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.spin = m
        self.Ms = Ms
        self.field = np.zeros(3 * mesh.n)

        if self.chi == 0.0:
            self.chi_inv = 0
        else:
            self.chi_inv = 1.0 / self.chi

    def compute_field(self, t=0):

        clib.compute_relaxation_field(self.spin,
                                      self.field,
                                      self.Ms,
                                      self.chi_inv,
                                      self.mesh.n)

        return self.field

    def compute_energy(self):
        # Not an energy term: this relaxes |m| towards Ms and contributes
        # nothing to the energy. The attributes follow the same convention as
        # the interaction classes
        self.energy = np.zeros(self.mesh.n)
        self.total_energy = 0.0

        return self.total_energy


class Laplace:

    """
        compute the laplace for given field.
    """

    def __init__(self, mesh):
        self.dx = mesh.dx * mesh.unit_length
        self.dy = mesh.dy * mesh.unit_length
        self.dz = mesh.dz * mesh.unit_length
        self.nx = mesh.nx
        self.ny = mesh.ny
        self.nz = mesh.nz
        self.field = np.zeros(3 * mesh.n)

    def compute_laplace_field(self, h, Ms):

        clib.compute_laplace_field(h, self.field,
                                   Ms,
                                   self.dx,
                                   self.dy,
                                   self.dz,
                                   self.nx,
                                   self.ny,
                                   self.nz)

        return self.field
