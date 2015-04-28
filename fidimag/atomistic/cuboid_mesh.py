"""
Implement a cuboid mesh.

TODO:
1. Implement cuboid mesh.
2. Write exchange field that uses said mesh (copy and modify existing exchange)
3. Compare with existing code.

Difference between cuboid and hexagonal mesh:
    How to compute coordinates (needed for DMI)
    How to find nearest neighbours

Let's see if we can get away with one class for both types of mesh or not.
Starting with cuboid case should enable us to make fastest progress.

"""
import numpy as np


class CuboidMesh(object):
    def __init__(Lx, Ly, Lz, nx, ny, nz):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.dx = float(Lx) / nx
        self.dy = float(Ly) / ny
        self.dz = float(Lz) / nz

        vertices = (nx + 1) * (ny + 1) * (nz + 1)
        self.coordinates = np.zeros((vertices, 3))
        self.init_coordinates()

    def init_coordinates(self):
        # do we really want x to be first dimension?
        for i in xrange(self.nx):
            for j in xrange(self.ny):
                for k in xrange(self.nz):
                    r = (i * self.dx + self.dx / 2.0,
                         j * self.dy + self.dy / 2.0,
                         k * self.dz + self.dz / 2.0)
                    index = 0 # function of i, j, k
                    self.coordinates[index] = r

    def neighbours(self):
        pass
        # this should happen on the C level actually
