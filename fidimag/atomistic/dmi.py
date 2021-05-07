import fidimag.extensions.clib as clib
from .energy import Energy
import numpy as np


class DMI(Energy):

    """

    This class provides the Dzyaloshinskii-Moriya Interaction (DMI) energy
    term, defined as

                  __  ->         ->      ->
         E =     \    D_ij   * ( S_i  X  S_j )
                 /__
                <i, j>
                i != j

    where D_ij is the Dzyaloshinskii vector, S_i and S_j are the total spin
    vectors at the i-th and j-th lattice sites, and <i, j> means counting the
    interaction between neighbouring spins only once. The Dzyaloshinskii vector
    magnitude can vary with space.

    Currently, there are implemented two kinds of DMI, that depend on the
    Dzyaloshinskii vector definition. Calling D_i the Dzyaloshinskii vector
    magnitude at the i-th site:

    BULK:
        ->         ^
        D_ij = D_i r_ij   with r_ij as the unit vector pointing from the i-th
                          towards the j-th site
                          (only working for square lattices)

    INTERFACIAL:

        ->           ^       ^
        D_ij = D_i ( r_ij X  z ) with r_ij as the unit vector pointing from the
                                 i-th towards the j-th site and z the unitary
                                 vector perpendicular to the magnetic material
                                 plane which, in theory, is in contact with a
                                 material with larger SOC (e.g. a metal).
                                 Thus, this DMI is defined in 2D, for
                                 interfaces or thin films.

    For further details about the DMI calculations, take a look
    to the C library documentation (dmi.c)


    OPTIONAL ARGUMENTS: -------------------------------------------------------

        dmi_type        :: 'bulk' or 'interfacial'
        name            :: Interaction name

    USAGE: --------------------------------------------------------------------

    If the DMI is homogeneous, it can be specified in a simulation object
    *Sim* as

            Sim.add(DMI(D, dmi_type='bulk'))

    where D is a float.

    Otherwise, it can be specified as any Fidimag scalar field, passing a
    function or an array. For example, a space dependent DMI that changes
    linearly in the x-direction can be defined as:

        def my_DMI(pos):
            D = 0.01 * meV
            return D * pos[0]

        # Add DMI to Simulation object
        Sim.add(DMI(my_DMI))

    """

    def __init__(self, D, name='DMI', dmi_type='bulk'):
        self.D = D

        self.name = name
        self.dmi_type = dmi_type

        self.jac = True

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(DMI, self).setup(mesh, spin, mu_s, mu_s_inv)

        if self.mesh_type == 'hexagonal':
            self.n_ngbs_dmi = 6
            # rdim = 3

        elif self.mesh_type == 'cuboid':
            self.n_ngbs_dmi = 4
            # rdim = 3

        # We will generate the Dzyaloshinskii vectors according
        # to the lattice, for the Interfacial DMI
        self.DMI_vector = self.compute_DMI_vectors(self.n_ngbs_dmi)

        if self.dmi_type == 'bulk':
            self._D = np.zeros(self.neighbours.shape)
            if isinstance(self.D, (int, float)):
                self._D[:, :] = self.D
            elif hasattr(self.D, '__call__'):
                n = self.mesh.n
                for i in range(n):
                    value = self.D(self.coordinates[i])
                    if len(value) == 6:
                        self._D[i, :] = value[:]
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
                                   self.mu_s_inv,
                                   self.energy,
                                   self._D,
                                   self.neighbours,
                                   self.n,
                                   self.n_ngbs
                                   )

        elif self.dmi_type == 'interfacial':

            clib.compute_dmi_field_interfacial(m,
                                               self.field,
                                               self.mu_s_inv,
                                               self.energy,
                                               self._D,
                                               self.neighbours,
                                               self.n,
                                               self.n_ngbs,
                                               self.n_ngbs_dmi,
                                               self.DMI_vector
                                               )

        return self.field

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

    def compute_DMI_vectors(self, n_ngbs_dmi):
        """

        Returns a DMI vectors array, according to the order of the nearest
        neighbours. So far we have a cubic and a hexagonal lattice, whose DMI
        are computed in 2D. The NNs arrays are in the order:

        Cubic       :: 4 NNs --> [-x, +x, -y, +y]
        Hexagonal   :: 6 NNs --> [right, left,
                                  top right, bottom left,
                                  top left, bottom right]

        Then, the DMI vectors array has (3 * n_ngbs_dmi) entries, 3 for every
        neighbouring site The vectors are normalised and computed according to:
        D_ij = r_ij X z where r_ij is the vector connecting a lattice site with
        the j-th neighbour (see Rohart and Thiaville PRB 88, 184422)

        """

        dmi_vec = np.zeros(n_ngbs_dmi * 3).reshape(-1, 3)
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
        for j in range(n_ngbs_dmi):
            dmi_vec[j] = np.cross(r_vec[j], z_vec)
            dmi_vec[j] /= np.linalg.norm(dmi_vec[j])

        return dmi_vec.flatten()
