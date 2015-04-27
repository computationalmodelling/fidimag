import numpy as np


def extract_mesh_info(dx, nx, x0, x1):
    if x1 is not None:
        nx = 1.0 * (x1 - x0) / dx
        if abs(nx - int(nx)) != 0:
            raise Exception('mesh parameters error!')
        return 1.0 * dx, int(nx), 1.0 * x0, 1.0 * x1
    else:
        x1 = 1.0 * dx * nx + x0
        return 1.0 * dx, int(nx), 1.0 * x0, x1


class FDMesh():

    def __init__(self, dx=1, dy=1, dz=1, nx=1, ny=1, nz=1, x0=0, y0=0, z0=0, x1=None, y1=None, z1=None, unit_length=1.0, pbc=None):
        """
        pbc could be None, 'x', 'y' or 'xy'.
        """

        self.dx, self.nx, self.x0, self.x1 = extract_mesh_info(dx, nx, x0, x1)
        self.dy, self.ny, self.y0, self.y1 = extract_mesh_info(dy, ny, y0, y1)
        self.dz, self.nz, self.z0, self.z1 = extract_mesh_info(dz, nz, z0, z1)

        self.dx_real = self.dx * unit_length
        self.dy_real = self.dy * unit_length
        self.dz_real = self.dz * unit_length

        self.nyz = self.ny * self.nz
        self.nxyz = self.nx * self.nyz
        self.unit_length = unit_length

        self.pbc = pbc

        self.compute_pos()

        self.cellsize = self.dx * self.dy * self.dz * unit_length**3

        self.xperiodic = 0
        self.yperiodic = 0

        if pbc is None:
            pass
        elif pbc == 'x' or pbc == '1d':
            self.xperiodic = 1
        elif pbc == 'y':
            self.yperiodic = 1
        elif pbc == 'xy' or pbc == '2d':
            self.xperiodic = 1
            self.yperiodic = 1
        else:
            raise Exception(
                "Only options None, 'x' ('1d'), 'y' or 'xy' ('2d') are acceptable!")

    def compute_pos(self):
        self.pos = []
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):

                    tp = ((i + 0.5) * self.dx + self.x0,
                          (j + 0.5) * self.dy + self.y0,
                          (k + 0.5) * self.dz + self.z0)

                    self.pos.append(tp)

    def index(self, i, j, k):
        idx = i * self.nyz + j * self.nz + k
        return idx

    def index_z(self, k=0):
        ids = [self.index(i, j, k) for i in range(self.nx)
               for j in range(self.ny)]
        return np.array(ids)

    # only used for tests
    def pos_at(self, i, j, k):
        idx = self.index(i, j, k)
        return self.pos[idx]
