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


def vector_from_coordinates(f, mesh):
    """
    Returns a np.array with values specified by `f`, where `f` should
    be a function that returns a vector when getting the coordinates of a
    cell of `mesh`.

    """
    v = np.zeros(mesh.vector_shape())
    for i, r in enumerate(mesh.coordinates):
        v[i] = f(r)
    return v


def scalar_from_coordinates(f, mesh):
    """
    Returns a np.array with values specified by `f`, where `f` should
    be a function that returns a number when getting the coordinates of a
    cell of `mesh`.

    """
    v = np.zeros(mesh.scalar_shape())
    for i, r in enumerate(mesh.coordinates):
        v[i] = f(r)
    return v


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
                    neighbours = set([cell for cell in [
                        self.index(i + 1, j, k),  # over
                        self.index(i - 1, j, k),  # under
                        self.index(i, j + 1, k),  # behind
                        self.index(i, j - 1, k),  # in front
                        self.index(i, j, k + 1),  # right
                        self.index(i, j, k - 1),  # left
                    ] if cell is not False])
                    connectivity.append(neighbours)
        return connectivity

    def index(self, i, j, k):
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
