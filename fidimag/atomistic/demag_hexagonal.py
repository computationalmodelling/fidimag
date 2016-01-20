import fidimag.extensions.dipolar as clib
import numpy as np
from fidimag.common import CuboidMesh


class DemagHexagonal(object):
    """
    This class allows to compute the Demag in a hexagonal mesh, using
    an equivalent cuboid mesh, where zeroes are appended between
    consecutive lattice sites in the original one

        The original mesh looks like

                6           7         8
                o -- -- -- o -- -- -- o

        3  o-- -- -- o -- -- -- o

                o -- -- -- o -- -- -- o
                0          1         2


        For every lattice site (--o-- in the figure), the scalar or vector
        field has associated a single value or a 3-vector.  Thus, we want to
        convert this mesh to a cuboid one, adding zero values between
        consecutive sites (spins) as follows:

            11   12                  16
            x -- o -- x -- o -- x -- o

         6  o -- x -- o -- x -- o -- x

            x -- o -- x -- o -- x -- o
            0    1    2    3    4   5

        Things to notice:
            * Every row has now the double number of nodes in the x direction,
              i.e. 2 * nx
            * The lattice spacing between cosequtive sites are half of those
              in the hexagonal lattice x (the y distance is the same),
              hence dx_c = 0.5 * dx
            * Even rows zeroes start from the 0th column while in odd
              rows, the start from the 1st column

        To convert scalar or vector fields to the new cuboid mesh, two
        functions are defined: scalar2cuboid and vector2cuboid.

    """

    def __init__(self, name='demag_hex'):
        self.name = name
        self.jac = True

    def setup(self, mesh, spin, mu_s):

        if mesh.mesh_type != 'hexagonal':
            raise Exception('This interaction is only defined'
                            'for hexagonal meshes'
                            )
        self.mesh = mesh

        self.dx = mesh.dx
        self.dx_c = 0.5 * mesh.dx

        self.dy = mesh.dy
        self.dz = mesh.dz

        self.nx = mesh.nx
        self.nx_c = 2 * mesh.nx

        self.ny = mesh.ny
        self.nz = mesh.nz

        self.spin = spin
        self.spin_c = np.zeros(2 * len(spin))

        self.n = mesh.n
        self.n_c = mesh.n * 2
        # self.n_c = self.mesh_c.n

        self.field = np.zeros(3 * self.n, dtype=np.float)
        self.field_c = np.zeros(3 * self.n_c, dtype=np.float)

        unit_length = mesh.unit_length
        self.mu_s_scale = np.zeros(mesh.n, dtype=np.float)
        self.mu_s_scale_c = np.zeros(2 * self.n, dtype=np.float)

        # note that the 1e-7 comes from \frac{\mu_0}{4\pi}
        self.scale = 1e-7 / unit_length ** 3

        # could be wrong, needs carefully tests!!!
        self.mu_s_scale = mu_s * self.scale

        # This seems not to be necessary
        # self.create_cuboid_mesh()

        self.scalar2cuboid(self.mu_s_scale, self.mu_s_scale_c)

        # Compute demag using the cuboid mesh with the double
        # number of x finite differences
        self.demag = clib.FFTDemag(self.dx_c, self.dy, self.dz,
                                   self.nx_c, self.ny, self.nz,
                                   tensor_type='dipolar')

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            # Copy the spin components to the cuboid mesh system
            self.vector2cuboid(spin, self.spin_c)
        else:
            # Copy the spin components to the cuboid mesh system
            self.vector2cuboid(self.spin, self.spin_c)

        m = self.spin_c

        # self.demag.compute_field(m, self.mu_s_scale, self.field)
        self.demag.compute_field(m, self.mu_s_scale_c, self.field_c)

        self.vector2cuboid(self.field, self.field_c, invert=True)

        # Return the hexagonal mesh field using the values
        # from the cuboid field
        return self.field

    def compute_exact(self):
        field = np.zeros(3 * self.n)
        field_c = np.zeros(2 * 3 * self.n)

        self.demag.compute_exact(self.spin_c, self.mu_s_scale_c, field_c)

        self.vector2cuboid(field, field_c, invert=True)

        return field

    def compute_energy(self):

        # Original:
        # energy = self.demag.compute_energy(
        #     self.spin, self.mu_s_scale, self.field)

        # We don't need to convert since energy is only a scalar
        self.vector2cuboid(self.field, self.field_c)
        self.vector2cuboid(self.spin, self.spin_c)

        energy = self.demag.compute_energy(
            self.spin_c, self.mu_s_scale_c, self.field_c)

        return energy / self.scale

    def create_cuboid_mesh(self):
        """

        This function will create a cuboid mesh based on the hexagonal one,
        with lattice points in between consecutive sites with mu_s = 0

        """
        nx, ny, nz = 2 * self.nx, self.ny, self.nz
        dx = self.dx * 0.5
        dy, dz = self.dy, self.dz

        mesh_c = CuboidMesh(nx, ny, nz,
                            dx, dy, dz,
                            unit_length=1e-9
                            )

        self.mesh_c = mesh_c

    def scalar2cuboid(self, scalar_field, scalar_field_c, invert=False):
        """

        This function copies a scalar field entries from a hexagonal mesh, into
        the equivalent cuboid mesh system.

        A scalar field array has the structure: f0, f1, f2, ...
        The cuboid scalar field has 2 times the length of the original
        scalar field array

        For this case, we:
            * Reshape both, the original and cuboid arrays with rows
              of nx and 2 * nx of length, respectively
            * For even rows [::2] , we copy the elements from the
              original mesh, to the cuboid field array, by
              adding the values every 2 columns, skipping the first one (0th),
              i.e. odd columns acquire the scalar field values, this is
              achieved with the indexing: [:, 1::2]
            * For odd rows [1::2], the process is similar, only that we do
              do not skip the 0th column, i.e. we do [:, ::2]

        """

        # Reshape both scalar field according to the number of elements
        # in the x direction
        scalar_field_c = scalar_field_c.reshape(-1, 2 * self.nx)
        scalar_field = scalar_field.reshape(-1, self.nx)

        if not invert:
            scalar_field_c[::2][:, 1::2] = scalar_field[::2]
            scalar_field_c[1::2][:, ::2] = scalar_field[1::2]
        else:
            scalar_field[::2] = scalar_field_c[::2][:, 1::2]
            scalar_field[1::2] = scalar_field_c[1::2][:, ::2]

        scalar_field_c = scalar_field_c.reshape(-1, )
        scalar_field = scalar_field.reshape(-1, )

    def vector2cuboid(self, vector_field, vector_field_c, invert=False):
        """

        This function copies a vector field entries from a hexagonal mesh, into
        the equivalent cuboid mesh system.

        A vector field array has the structure:
                fx0, fy0, fz0, fx1, fy1, fz1, fx2, ...

        The cuboid vector field has 2 times the length of the original
        vector field array.

        The procedure is similar to the one in the scalar field, the only
        difference is in the reshaping of the vector field arrays. We
        make it such that instead of a single value, we use 3 vectors, thus
        we have a 3 dimensional matrix where every row has still nx or 2 * nx
        values.

        For example, in a 2 * 3 system:

        Original: vector_field = [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                                   11, 12, 13, 14, 15, 16, 17]

        Reshaped: [ [ [ 0,  1,  2], [ 3,  4,  5] ],

                    [ [ 6,  7,  8], [ 9, 10, 11] ],

                    [ [12, 13, 14], [15, 16, 17] ]
                  ]

                  So we have 2 sites per row,as the original system, where
                  a 3-vector is associated to each entry, making a
                  (2, 3, 3) matrix --> (ny, nx, 3) in general

                  Thus, for the cuboid version, we would have the double
                  of columns, i.e. a (ny, 2 * nx, 3) matrix and copy
                  the elements as in the scalar field

        """
        vector_field_c = vector_field_c.reshape(-1, 2 * self.nx, 3)
        vector_field = vector_field.reshape(-1, self.nx, 3)

        if not invert:
            vector_field_c[::2][:, 1::2] = vector_field[::2]
            vector_field_c[1::2][:, ::2] = vector_field[1::2]
        else:
            vector_field[::2] = vector_field_c[::2][:, 1::2]
            vector_field[1::2] = vector_field_c[1::2][:, ::2]

        vector_field_c = vector_field_c.reshape(-1,)
        vector_field = vector_field.reshape(-1)
