import fidimag.extensions.clib as clib
from energy import Energy
import numpy as np


class DMI(Energy):

    """
    Hamiltonian = D*[S_i x S_j]
        ==> H = D x S_j
    """

    def __init__(self, D, name='dmi', dmi_type='bulk'):
        self.D = D

        self.name = name
        self.dmi_type = dmi_type

        self.jac = True

    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        if self.dmi_type == 'bulk':
            clib.compute_dmi_field(m,
                                   self.field,
                                   self.energy,
                                   self.D,
                                   self.neighbours,
                                   self.n)

        elif self.dmi_type == 'interfacial':

            # We will generate the Dzyaloshinskii vectors according
            # to the lattice, for the Interfacial DMI
            if self.mesh_type == 'hexagonal':
                self.nneighbours = 6
                # rdim = 3

            elif self.mesh_type == 'cuboid':
                self.nneighbours = 4
                # rdim = 3

            self.DMI_vector = self.compute_DMI_vectors(self.nneighbours)

            clib.compute_dmi_field_interfacial(m,
                                               self.field,
                                               self.energy,
                                               self.D,
                                               self.neighbours,
                                               self.n,
                                               self.nneighbours,
                                               self.DMI_vector
                                               )

        return self.field * self.mu_s_inv

    def compute_energy_direct(self):
        """
        mainly for testing
        """
        energy = clib.compute_dmi_energy(self.spin,
                                         self.D,
                                         self.nx,
                                         self.ny,
                                         self.nz,
                                         self.xperiodic,
                                         self.yperiodic)

        return energy

    def compute_DMI_vectors(self, nneighbours):

        dmi_vec = np.zeros(nneighbours * 3).reshape(-1, 3)
        r_vec = np.zeros_like(dmi_vec)
        z_vec = np.array([0, 0, 1])

        rij = lambda i, j: np.array([i * self.mesh.dx * 0.5,
                                     j * self.mesh.dy * 3.0 / 4.0,
                                     0])

        if self.mesh_type == 'hexagonal':
            r_vec = np.array([rij(2, 0), -rij(2, 0),   # right and left
                              rij(1, 1), -rij(1, 1),   # top right and btm left
                              rij(-1, 1), -rij(-1, 1)  # top left and btm right
                              ])

        elif self.mesh_type == 'cuboid':
            r_vec = np.array([[-1, 0, 0], [1, 0, 0],
                              [0, -1, 0], [0, 1, 0]
                              ])

        for j in range(nneighbours):
            dmi_vec[j] = np.cross(r_vec[j], z_vec)
            dmi_vec[j] /= np.linalg.norm(dmi_vec[j])

        return dmi_vec.flatten()
