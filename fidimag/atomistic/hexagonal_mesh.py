
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
        self.nz = 1  # time will tell if 0 is a better value here
        self.periodicity = periodicity

        self.dx = sqrt(3) * radius
        self.dy = 2.0 * radius

        # To avoid moodifying the other classes that assume a 3D sample
        self.dz = 1

        self.Lx = self.nx * self.dx
        self.Ly = self.ny * self.dy * 3.0 / 4.0 + self.dy / 4.0

        self.n = nx * ny  # total number of cells

        self.size = (self.nx, self.ny, 1)  # number of cells in all directions
        self.coordinates = self.init_coordinates()
        self.neighbours = self.init_neighbours()
        self.vertices, self.hexagons = self.init_grid()
        self.mesh_type = 'hexagonal'
        self.unit_length = unit_length

    def init_coordinates(self):
        coordinates = np.zeros((self.n, 3))
        for j in xrange(self.ny):
            for i in xrange(self.nx):
                index = self.index(i, j)
                r = (j * self.dx / 2.0 + i * self.dx + self.dx / 2.0,
                     j * self.dy * 3.0 / 4.0 + self.dy / 2.0,
                     0
                     )
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
                ]]
                neighbours = [other if other != cell
                              else -1 for other in neighbours]
                connectivity.append(neighbours)
        return np.array(connectivity, dtype=np.int32)

    def init_grid(self):
        """
        Compute the coordinates of the vertices that make up the hexagonal
        grid and how the vertices must be assembled to form hexagons.

        """
        vertex_counter = 0
        vertices = []  # list of tuples of coordinates
        hexagons = []  # list of tuples of vertices
        for j in xrange(self.ny):
            for i in xrange(self.nx):
                index = self._index(i, j)
                x, y = self.coordinates[index]
                corners = hexagon_corners(x, y, self.radius)
                hexagon = []
                # We'll go through the corners in a counter-clockwise direction.
                # For each corner, we think about if it's a "new" vertex, or
                # if it has been created by a neighbouring hexagon before.
                # NE (0) and N (1) corners will always be "new" vertices
                for c in (0, 1):
                    vertices.append(corners[c])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # NW (2) corner could have been created by western neighbour
                # where it will have been the the NE (0) corner
                W = self._index(i - 1, j)
                if W is not False:  # can't replace by if W because 0 == False
                    hexagon.append(hexagons[W][0])  # our NW (2) is west's NE (0)
                else:
                    vertices.append(corners[2])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # SW (3) corner could have been created either by western
                # or south-western neighbour
                SW = self._index(i, j - 1)
                if W is not False:
                    hexagon.append(hexagons[W][5])  # our SW is west's SE (5)
                elif SW:
                    hexagon.append(hexagons[SW][1])  # or south-west's N (1)
                else:
                    vertices.append(corners[3])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # S (4) corner could have been created by south-western neighbour
                if SW is not False:
                    hexagon.append(hexagons[SW][0])  # our S is south-west's NE (0)
                else:
                    vertices.append(corners[4])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # SE (5) corner could have been created by south-eastern neighbour
                SE = self._index(i + 1, j - 1)
                if SE is not False:
                    hexagon.append(hexagons[SE][1])  # our SE is south-east's N (1)
                else:
                    vertices.append(corners[5])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                hexagons.append(hexagon)
        return np.array(vertices), np.array(hexagons)

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
            return -1
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
