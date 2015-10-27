
"""
Represent hexagonal 2D mesh.

The hexagons are pointy topped (as against flat topped). We use axial
coordinates, also known as trapezoidal coordinates (compared to cubic or
offset coordinates).

"""
import numpy as np
from math import sqrt


class HexagonalMesh(object):
    def __init__(self, radius, nx, ny,
                 periodicity=(False, False), unit_length=1.0):
        """
        Create mesh with nx cells in x-direction and ny cells in the
        y-direction. The size of a hexagon is given by

            width  = sqrt(3) * radius
            height = 2 * radius.

        By default, the mesh is not periodic along any axis, so the periodicity
        is set to (False, False). Passing a tuple with any combination
        of entries set to True will enable periodicity along the given axes.

        """
        self.nx = nx
        self.ny = ny
        self.periodicity = periodicity

        self.dx = sqrt(3) * radius
        self.dy = 2.0 * radius

        self.Lx = self.nx * self.dx
        self.Ly = self.ny * self.dy * 3.0 / 4.0 + self.dy / 4.0

        self.n = nx * ny  # total number of cells

        self.coordinates = self.init_coordinates()
        self.neighbours = self.init_neighbours()

        self.mesh_type = 'hexagonal'
        self.unit_length = unit_length
        self.nxyz = 0

    def init_coordinates(self):
        coordinates = np.zeros((self.n, 2))
        for j in xrange(self.ny):
            for i in xrange(self.nx):
                index = self.index(i, j)
                r = (j * self.dx / 2.0 + i * self.dx + self.dx / 2.0,
                     j * self.dy * 3.0 / 4.0 + self.dy / 2.0)
                coordinates[index] = r
        return coordinates

    def init_neighbours(self):
        # can't use numpy array because not all cells have the same
        # number of neighbours (the cells on the boundaries have less)
        connectivity = []
        for j in xrange(self.ny):
            for i in xrange(self.nx):
                cell = self._index(i, j)
                # there can be periodicity along x and y axes, but not
                # along the third cubic axis (would be z), hence, we use
                # _index to disregard periodicity in the NW and SE directions.
                neighbours = [other for other in [
                    self.index(i + 1, j),       # right  east
                    self.index(i - 1, j),       # left   west
                    self.index(i, j + 1),       # up     north-east
                    self.index(i, j - 1),       # down   south-west
                    self._index(i - 1, j + 1),  #        north-west
                    self._index(i + 1, j - 1),  #        south-east
                ] if other is not False and other != cell]
                connectivity.append(neighbours)
        return connectivity

    def index(self, i, j):
        """
        Returns the index for the cell with ordinals i, j or False if that
        cell would be out of bounds. Handles periodic meshes.

        """
        if self.periodicity[0]:  # if mesh is periodic in x-direction
            if i == -1:          # then wrap the left side
                i = self.nx - 1  # to the right
            if i == self.nx:     # and wrap the right side
                i = 0            # to the left
        if self.periodicity[1]:
            if j == -1:
                j = self.ny - 1
            if j == self.ny:
                j = 0
        return self._index(i, j)

    def _index(self, i, j):
        """
        Returns the index for the cell with ordinals i, j
        or False if that cell would be out of bounds.

        """
        if i < 0 or j < 0 or i >= self.nx or j >= self.ny:
            return False
        return i + j * self.nx

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
