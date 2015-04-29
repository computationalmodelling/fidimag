import numpy as np


class CuboidMesh(object):
    def __init__(self, Lx, Ly, Lz, nx, ny, nz):
        """
        Create mesh with dimensions Lx, Ly, Lz.

        Divided into nx cells along the x-axis, ny cells along the y-axis
        and nz cells along the z-axis for a total of nx * ny * nz cells.

        """
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.dx = float(Lx) / nx  # size of one cell
        self.dy = float(Ly) / ny
        self.dz = float(Lz) / nz

        self.n = nx * ny * nz  # total number of cells
        self.nxy = nx * ny  # number of cells in the x-y plane

        self.cells = self.n
        self.coordinates = np.zeros((self.cells, 3))
        self.init_coordinates()

    def init_coordinates(self):
        for i in xrange(self.nz):
            for j in xrange(self.ny):
                for k in xrange(self.nx):
                    index = i * self.nxy + j * self.nx + k
                    r = (k * self.dx + self.dx / 2.0,
                         j * self.dy + self.dy / 2.0,
                         i * self.dz + self.dz / 2.0)
                    self.coordinates[index] = r
