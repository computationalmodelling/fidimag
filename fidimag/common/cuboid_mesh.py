"""
Represent cuboid mesh.

Indexing follows the x-axis first and then the y-axis before moving
up to the next horizontal slice, as shown in the sketch below. We have
chosen this indexing scheme because systems like thin films are usually
longest along the x-axis.

        +-------+
     .'       .:|
    +-------+:::|
    |       |:::|
    |   30  |::;+-------+-------+-------+-------+-------+
    |       |;'       .:| 11  .'  12  .'  13  .'  14  .:|
    +-------+-------+:::|---+-------+-------+-------+:::|
    |       |       |:::| .'   7  .'   8  .'   9  .:|:::|
    |   15  |  16   |::;+-------+-------+-------+:::|:::+
    |       |       |;'       .'      .'      .:|:::|::'
    +-------+-------+-------+-------+-------+:::|:::+'
    |       |       |       |       |       |:::|:.'
    |   0   |   1   |   2   |   3   |   4   |:::+'
    |       |       |       |       |       |::'
    +-------+-------+-------+-------+-------+'

N.B. This means that when iterating over the cells, the x-axis is traversed
in the innermost loop, and the z-axis in the outermost loop!

"""
import numpy as np


class CuboidMesh(object):
    def __init__(self, dx=1, dy=1, dz=1, nx=1, ny=1, nz=1,
                 periodicity=(False, False, False), unit_length=1.0):
        """
        Create mesh with cells of size dx * dy * dz.

        There will be nx cells along the x-axis, ny cells along the y-axis
        and nz cells along the z-axis for a total of nx * ny * nz cells.

        By default, the mesh is not periodic along any axis, so the periodicity
        is set to (False, False, False). Passing a tuple with any combination
        of entries set to True will enable periodicity along the given axes.

        Usage:
            mesh = CuboidMesh(2, 2, 2, 250, 25, 2, periodicity=(True, False, False))
            # create a mesh of dimensions 500 x 50 x 4 nm, with cellsize
            # of 2 nm in any direction and periodic along the x-axis.

        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.periodicity = periodicity

        self.Lx = dx * nx  # total size of mesh
        self.Ly = dy * ny
        self.Lz = dz * nz

        self.n = nx * ny * nz  # total number of cells
        self.nxy = nx * ny  # number of cells in the x-y plane

        self.mesh_type = "cuboid"
        self.unit_length = unit_length

        self.coordinates = self.init_coordinates()
        self.neighbours = self.init_neighbours()

    def init_coordinates(self):
        coordinates = np.zeros((self.n, 3))
        for i in xrange(self.nz):
            for j in xrange(self.ny):
                for k in xrange(self.nx):
                    index = self.index(i, j, k)
                    r = (k * self.dx + self.dx / 2.0,
                         j * self.dy + self.dy / 2.0,
                         i * self.dz + self.dz / 2.0)
                    coordinates[index] = r
        return coordinates

    def init_neighbours(self):
        # can't use numpy array because not all cells have the same
        # number of neighbours (the cells on the boundaries have less)
        connectivity = []
        for i in xrange(self.nz):
            for j in xrange(self.ny):
                for k in xrange(self.nx):
                    cell = self._index(i, j, k)
                    neighbours = [other for other in [
                        self.index(i + 1, j, k),  # over
                        self.index(i - 1, j, k),  # under
                        self.index(i, j + 1, k),  # behind
                        self.index(i, j - 1, k),  # in front
                        self.index(i, j, k + 1),  # right
                        self.index(i, j, k - 1),  # left
                    ] if other is not False and other != cell]
                    connectivity.append(neighbours)
        return connectivity

    def index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds. Handles periodic meshes.

        """
        if self.periodicity[0]:  # if mesh is periodic in x-direction
            if k == -1:          # then wrap the left side
                k = self.nx - 1  # to the right
            if k == self.nx:     # and wrap the right side
                k = 0            # to the left
        if self.periodicity[1]:
            if j == -1:
                j = self.ny - 1
            if j == self.ny:
                j = 0
        if self.periodicity[2]:
            if i == -1:
                i = self.nz - 1
            if i == self.nz:
                i = 0
        return self._index(i, j, k)

    def _index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds.

        """
        if i < 0 or j < 0 or k < 0 or i >= self.nz or j >= self.ny or k >= self.nx:
            return False
        return i * self.nxy + j * self.nx + k

    def cells(self):
        """
        Generator to iterate over the cell indices.

        """
        current = 0
        while current < self.n:
            yield current
            current += 1

    def scalar_shape(self):
        """
        Return the appropriate shape for np.array for scalar field over the cells.

        Usage example:
            alpha = np.zeros(mesh.scalar_shape())

        """
        return (self.n, 1)

    def vector_shape(self):
        """
        Return the appropriate shape for np.array for vector field over the cells.

        Usage example:
            m = np.zeros(mesh.vector_shape())

        """
        return (self.n, 3)
