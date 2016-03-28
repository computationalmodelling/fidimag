import fidimag.extensions.clib as clib
from .energy import Energy
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

    def setup(self, mesh, spin, mu_s):
        super(DMI, self).setup(mesh, spin, mu_s)

        if self.mesh_type == 'hexagonal':
            self.nneighbours = 6
            # rdim = 3

        elif self.mesh_type == 'cuboid':
            self.nneighbours = 4
            # rdim = 3

        # We will generate the Dzyaloshinskii vectors according
        # to the lattice, for the Interfacial DMI
        self.DMI_vector = self.compute_DMI_vectors(self.nneighbours)
        
        if self.dmi_type == 'bulk':
            self._D = np.zeros(self.neighbours.shape)
            if isinstance(self.D, (int, float)):
                self._D[:,:] = self.D
            elif hasattr(self.D, '__call__'):
                n =  self.mesh.n
                for i in range(n):
                    value = self.D(self.coordinates[i])
                    if len(value) == 6:
                        self._D[i,:] = value[:]
                    else:
                        raise Exception('The given spatial function for D is not acceptable!')
            

    def compute_field(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        if self.dmi_type == 'bulk':
              clib.compute_dmi_field(m,
                                   self.field,
                                   self.energy,
                                   self._D,
                                   self.neighbours,
                                   self.n)

        elif self.dmi_type == 'interfacial':

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
        """

        Returns a DMI vectors array, according to the order of the nearest
        neighbours. So far we have a cubic and a hexagonal lattice, whose DMI
        are computed in 2D. The NNs arrays are in the order:

        Cubic       :: 4 NNs --> [-x, +x, -y, +y]
        Hexagonal   :: 6 NNs --> [right, left,
                                  top right, bottom left,
                                  top left, bottom right]

        Then, the DMI vectors array has (3 * nneighbours) entries, 3 for every
        neighbouring site The vectors are normalised and computed according to:
        D_ij = r_ij X z where r_ij is the vector connecting a lattice site with
        the j-th neighbour (see Rohart and Thiaville PRB 88, 184422)

        """

        dmi_vec = np.zeros(nneighbours * 3).reshape(-1, 3)
        r_vec = np.zeros_like(dmi_vec)
        z_vec = np.array([0, 0, 1])

        rij = lambda i, j: np.array([i * self.mesh.dx * 0.5,
                                     j * self.mesh.dy,
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

        # Here we compute the cross product with the unitary z vector
        # and normalise
        for j in range(nneighbours):
            dmi_vec[j] = np.cross(r_vec[j], z_vec)
            dmi_vec[j] /= np.linalg.norm(dmi_vec[j])

        return dmi_vec.flatten()
