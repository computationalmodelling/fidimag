import fidimag.extensions.dipolar as clib
import numpy as np
from fidimag.common import CuboidMesh


class DemagHexagonal(object):

    def __init__(self, name='demag'):
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
        self.n = mesh.n
        self.n_c = self.mesh_c.n

        self.field = np.zeros(3 * self.n, dtype=np.float)
        self.field_c = np.zeros(3 * self.n_c, dtype=np.float)

        unit_length = mesh.unit_length
        self.mu_s_scale = np.zeros(mesh.n, dtype=np.float)

        # note that the 1e-7 comes from \frac{\mu_0}{4\pi}
        self.scale = 1e-7 / unit_length ** 3

        # could be wrong, needs carefully tests!!!
        self.mu_s_scale = mu_s * self.scale

        self.create_cuboid_mesh()

        self.demag = clib.FFTDemag(self.dx_c, self.dy, self.dz,
                                   self.nx_c, self.ny, self.nz,
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

    def create_cuboid_mesh(self):
        """

        This function will create a cuboid mesh based on the hexagonal one,
        with lattice points in between consequtive sites with mu_s = 0

        """
        nx, ny, nz = 2 * self.nx, self.ny, self.nz
        dx = self.dx * 0.5
        dy, dz = self.dy, self.dz

        mesh_c = CuboidMesh(nx, ny, nz,
                            dx, dy, dz,
                            unit_length=1e-9
                            )

        self.mesh_c = mesh_c

        self.mu_s_scale_c = np.zeros(2 * self.n, dtype=np.float)
        self.mu_s_scale_c = self.mu_s_scale_c.reshape(-1, 2 * self.nx)

        self.mu_s_scale = self.mu_s_scale.reshape(-1, self.nx)

        self.mu_s_scale_c[::2][:, ::2] = self.mu_s_scale[::2]
        self.mu_s_scale_c[1::2][:, 1::2] = self.mu_s_scale[1::2]

        self.mu_s_scale_c = self.mu_s_scale_c.reshape(-1, )
        self.mu_s_scale = self.mu_s_scale.reshape(-1, )

    def scalar2cuboid(self, scalar_field, scalar_field_c):

        scalar_field_c = self.mu_s_scale_c.reshape(-1, 2 * self.nx)

        self.mu_s_scale = self.mu_s_scale.reshape(-1, self.nx)

        self.mu_s_scale_c[::2][:, ::2] = self.mu_s_scale[::2]
        self.mu_s_scale_c[1::2][:, 1::2] = self.mu_s_scale[1::2]

        self.mu_s_scale_c = self.mu_s_scale_c.reshape(-1, )
        self.mu_s_scale = self.mu_s_scale.reshape(-1, )
