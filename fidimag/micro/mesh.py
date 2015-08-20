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

        self.nxy = self.nx*self.ny
        self.nxyz = self.nxy * self.nz
        self.unit_length = unit_length

        self.pbc = pbc

        self.compute_pos()

        self.cellsize = self.dx * self.dy * self.dz * unit_length**3

        self.xperiodic = 0
        self.yperiodic = 0
        self.zperiodic = 0

        if pbc is None:
            pass
        elif pbc == '1d':
            self.xperiodic = 1
        elif pbc == '2d':
            self.xperiodic = 1
            self.yperiodic = 1
        else:
            if 'x' in pbc:
                self.xperiodic = 1
            if 'y' in pbc:
                self.yperiodic = 1
            if 'z' in pbc:
                self.zperiodic = 1

        self.init_neighbours()

    def compute_pos(self):

        self.pos = []
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):

                    tp = ((i + 0.5) * self.dx + self.x0,
                          (j + 0.5) * self.dy + self.y0,
                          (k + 0.5) * self.dz + self.z0)

                    self.pos.append(tp)

    def index(self, i, j, k):
        idx = k * self.nxy + j * self.nx + i
        return idx

    def init_neighbours(self):
        """

        Creates a *connectivity* array with the index of the neighbours for
        every lattice site in the order:

               | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  |
                 i=0                           i=1                ...

        where  0-y  is the index of the neighbour of the 0th spin,
        in the -y direction, for example. So we can access the -x neighbour
        of the ith spin using 6*i, the +x neighbour with 6*i+1,
        the -y neighbour with 6*i+2, and so on

        """
        connectivity = []

        for k in xrange(self.nz):
            for j in xrange(self.ny):
                for i in xrange(self.nx):
                    ngbs = [
                        self.index(i - 1, j, k),  # x-1
                        self.index(i + 1, j, k),  # x+1
                        self.index(i, j - 1, k),  # y-1
                        self.index(i, j + 1, k),  # y+1
                        self.index(i, j, k - 1),  # z-1
                        self.index(i, j, k + 1),  # z+1
                    ]
                    connectivity.append(ngbs)

        self.connectivity = np.array(connectivity, dtype=np.int32)

    def index_z(self, k=0):
        ids = [self.index(i, j, k) for i in range(self.nx)
               for j in range(self.ny)]
        return np.array(ids)

    # only used for tests
    def pos_at(self, i, j, k):
        idx = self.index(i, j, k)
        return self.pos[idx]
