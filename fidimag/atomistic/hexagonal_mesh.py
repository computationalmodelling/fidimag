
"""
Represent hexagonal 2D mesh.

The hexagons are pointy topped (as against flat topped). We use axial
coordinates, also known as trapezoidal coordinates (compared to cubic or
offset coordinates).

We have two different alignments of the hexagons:

    DIAGONAL                       and SQUARE

         /\ /\ /\                           /\ /\ /\
        |  |  |  |                         |  |  |  |
        | 6| 7| 8|                         | 6| 7| 8|        j = 2
       /\ /\ /\ /                         /\ /\ /\ /
      |  |  |  |                         |  |  |  |
      | 3| 4| 5|                         | 3| 4| 5|          j = 1
     /\ /\ /\ /                           \ /\ /\ /\
    |  |  |  |                             |  |  |  |
    | 0| 1| 2|                             | 0| 1| 2|        j = 0
     \/ \/ \/                               \/ \/ \/

For every case, it will be necessary to generate the neighbours matrix
using different index orders. In particular, for the square
arrangement, the ordering changes for even and odd numbered rows

In both cases, every row of the neighbours matrix has the same order
for the indexes:

    | left right top_right bottom_left top_left bottom_right |

Periodic boundary conditions are well defined in the x direction, however,
along the y axis, it is only well defined for the DIAGONAL case, where,
for example, we have periodicity between 0 and 6, 1 and 7 and 2 and 8.
For the SQUARE case, it is still not properly defined.

Coordinates a

"""
import numpy as np
from math import sqrt

from six.moves import range


class HexagonalMesh(object):
    def __init__(self, radius, nx, ny,
                 periodicity=(False, False),
                 unit_length=1.0,
                 alignment='diagonal',
                 shells=1):
        """
        Create mesh with nx cells in x-direction and ny cells in the
        y-direction. The size of a hexagon is given by

            height  = sqrt(3) * radius
            width = 2 * radius.

        radius ::

        The radius is for an incircle inside the hexagon, thus the distance
        between two lattice sites at the same y coordinate, is just 2 * radius,
        which can be seen as the triangular lattice constant.
        The height refers to the y distance between two lattice sites of
        consecutive rows.
        The first row of atoms has the lattice points located at the middle
        of an hexagon, whose height is not the same distance than
        the *height* defined before. This
        distance, in terms of dx = 2 * radius, can be easily calculated as
        2 * dx / sqrt(3), since the hexagon height is: (3 / 4) * height

        alignment ::

        The alignment of the hexagons can be set to 'diagonal' or 'square'. In
        both cases the matrix with the neighbours indexes will have the same
        order for every row (see definition of shells)

        periodicity ::

        By default, the mesh is not periodic along any axis, so the periodicity
        is set to (False, False). Passing a tuple with any combination of
        entries set to True will enable periodicity along the given axes.
        Periodic boundaries are not defined when using a 'square' alignment.

        unit_length ::

        The scale unit for the lattice distances. Default is 1.0 but it is
        commonly used a nm = 1e-9

        shells ::

        An integer specifying the number of shells of neighbours to be computed
        and stored in the *ngbs array. By default, a value of 1 indicates
        only nearest neighbours. The *ngbs array has the structure:

        [ [ ngbs_1st_shell_0 negbs_2nd_shell_0 ...  ],      --> ngbs of spin 0
          [ ngbs_1st_shell_1 negbs_2nd_shell_1 ...  ],      --> ngbs of spin 1
          ...
          [ ngbs_1st_shell_n negbs_2nd_shell_n ...  ]       --> ngbs of spin n
        ]

        where ngbs_1st_shell_i are the nearest neighbours of the i-th spin in
        the order:
            [ left right top_right bottom_left top_left bottom_right ]

        For the other shells see the *_ngbs_Xth_shell methods from this class

        """

        if shells > 9:
            raise ValueError('Number of shells cannot be larger than 8')
        else:
            self.n_shells = shells

        # Number of neighbours (ngbs) according to the shell, i.e. at the
        # first shell (nearest ngbs) there are 6 ngbs,
        # second shell (NNNs) -> 6 ngbs, etc
        # (we set zero to 0 ngbs to make for loops more understandable)
        self._n_ngbs_shell = np.array([0, 6, 6, 6, 12, 6, 6, 12, 6, 12],
                                      dtype=np.int32)

        # Total number of ngbs:
        self.n_ngbs = np.sum([self._n_ngbs_shell[i]
                              for i in range(1, self.n_shells + 1)])

        # List with the sum of number of neighbours, to set the range of cols
        # to store the ngbs indexes for a specific shell in a specific row. For
        # example, 1st ngbs are stored in cols 0-5, 2nd ngbs in 6-11, etc.
        self._sum_ngbs_shell = np.array([np.sum([self._n_ngbs_shell[i]
                                                 for i in range(max_sh + 1)])
                                         for max_sh in range(self.n_shells + 1)],
                                        dtype=np.int32)

        # Dictionary to call the methods that return the indexes of the
        # neighbours for a specific shell (like a switch statement)
        self._ngbs_i_shell = {1: self._ngbs_first_shell,
                              2: self._ngbs_second_shell,
                              3: self._ngbs_third_shell,
                              4: self._ngbs_fourth_shell,
                              5: self._ngbs_fifth_shell,
                              6: self._ngbs_sixth_shell,
                              7: self._ngbs_seventh_shell,
                              8: self._ngbs_eigth_shell,
                              9: self._ngbs_ninth_shell
                              }

        self.nx = nx
        self.ny = ny
        self.nz = 1  # time will tell if 0 is a better value here
        self.periodicity = periodicity

        self.dy = sqrt(3) * radius
        self.dx = 2.0 * radius
        self.radius = radius

        # To avoid moodifying the other classes that assume a 3D sample
        self.dz = 1

        # hexagons height: h = (3 / 4) * dy
        self.h = self.dx * 2. / np.sqrt(3)

        self.Lx = self.nx * self.dx
        # This is: (n - 1) * self.dy + self.h
        self.Ly = self.ny * self.dy + self.dy / 3.

        self.n = nx * ny  # total number of cells

        self.alignment = alignment

        self.size = (self.nx, self.ny, 1)  # number of cells in all directions
        self.coordinates = self.init_coordinates()
        self.neighbours = self.init_neighbours()
        # self.vertices, self.hexagons = self.init_grid()
        self.mesh_type = 'hexagonal'
        self.unit_length = unit_length

        self.vertices, self.hexagons = self.init_grid()

    def init_coordinates(self):
        coordinates = np.zeros((self.n, 3))
        for j in range(self.ny):
            for i in range(self.nx):
                index = self.index(i, j)
                # For a diagonal alignment, the hexagons are shifted
                # in their x-position by dx * 0.5, on every row
                if self.alignment == 'diagonal':
                    r = (j * self.dx / 2.0 + i * self.dx + self.dx / 2.0,
                         j * self.dy + self.h / 2.0,
                         0
                         )
                # For a square alignment, the hexagons will
                # always be in the same x-position for even numbered rows (j)
                # x(i=0) = dx * 0.5
                elif self.alignment == 'square':
                    if j % 2 == 0:
                        sign = 1
                    else:
                        sign = 0

                    r = (sign * self.dx / 2.0 + i * self.dx + self.dx / 2.0,
                         j * self.dy + self.h / 2.0,
                         0
                         )

                coordinates[index] = r
        return coordinates

    def init_neighbours(self):

        # The neighbours array will be a NxM array where N (rows) will be the
        # numer of lattice sites and M is the total number of ngbs per site
        neighbours = np.zeros((self.nx * self.ny, self.n_ngbs), dtype=np.int32)

        site = 0  # Counter for lattice sites so we populate array row wise

        # Sweep through every lattice site, thus every row of the neighbours
        # array
        # For every row, we will compute the neighbours indexes for the *sh*
        # shell, calling the corresponding function from the _ngbs_i_shell
        # dictionary.  We store the neighbours according to the shells order,
        # i.e. for the k-th site:
        #
        #                   1st shell ngbs  | 2nd shell ngbs  | 3rd shell ...
        # neighbours[k] = [x  x  x  x  x  x  o  o  o  o  o  o  -  -  -  -  -  - ...]
        #
        for j in range(self.ny):
            for i in range(self.nx):
                cell = self._index(i, j)

                for sh in range(1, self.n_shells + 1):
                    # Save ngbs indexes in the corresponding range of the row
                    ngbs_range = slice(self._sum_ngbs_shell[sh - 1],
                                       self._sum_ngbs_shell[sh])
                    # compute and store them in this range
                    neighbours[site, ngbs_range] = self._ngbs_i_shell[sh](i, j)

                    # Set indexes to -1 for sites outside the lattice
                    neighbours[site, ngbs_range] = [other if other != cell
                                                    else -1 for other in
                                                    neighbours[site, ngbs_range]]

                    # we must break to populate *neighbours* only according to
                    # the specified number of shells
                    if self.n_shells == sh:
                        continue

                site += 1

        return neighbours

    def init_grid(self):
        """
        Compute the coordinates of the vertices that make up the hexagonal
        grid and how the vertices must be assembled to form hexagons.

        """
        vertex_counter = 0
        vertices = []  # list of tuples of coordinates
        hexagons = []  # list of tuples of vertices
        for j in range(self.ny):
            for i in range(self.nx):
                index = self._index(i, j)
                x, y = self.coordinates[index][0], self.coordinates[index][1]
                # self.radius is the inradius while self.h/2  is the circumradius
                corners = self.hexagon_corners(x, y, self.h * 0.5)
                hexagon = []
                # We'll go through the corners in a counter-clockwise direction.
                # For each corner, we think about if it's a "new" vertex, or
                # if it has been created by a neighbouring hexagon before.

                # Here we define the neighbours indexes to check if vertexes
                # were already created
                if self.alignment == 'square':
                    if j % 2 == 0:
                        W = self.index(i - 1, j)        # left   west
                        SW = self.index(i, j - 1)       # down   south-west
                        SE = self.index(i + 1, j - 1)   #        south-east
                    else:
                        W = self.index(i - 1, j)        # left   west
                        SW = self.index(i - 1, j - 1)   # down   south-west
                        SE = self.index(i, j - 1)       #        south-east

                elif self.alignment == 'diagonal':
                    W = self._index(i - 1, j)
                    SW = self._index(i, j - 1)
                    SE = self._index(i + 1, j - 1)


                # NE (0) and N (1) corners will always be "new" vertices
                for c in (0, 1):
                    vertices.append(corners[c])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # NW (2) corner could have been created by western neighbour
                # where it will have been the the NE (0) corner
                # Sites with no magnetisation have a value of -1 (before it was
                # False but we changed to numpy arrays; I will
                # let the False statements just in case)
                if W is not (False or -1):  # can't replace by if W because 0 == False
                    hexagon.append(hexagons[W][0])  # our NW (2) is west's NE (0)
                else:
                    vertices.append(corners[2])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # SW (3) corner could have been created either by western
                # or south-western neighbour
                if W is not (False or -1):
                    hexagon.append(hexagons[W][5])  # our SW is west's SE (5)
                elif SW is not (False or -1):
                    hexagon.append(hexagons[SW][1])  # or south-west's N (1)
                else:
                    vertices.append(corners[3])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # S (4) corner could have been created by south-western neighbour
                if SW is not (False or -1):
                    hexagon.append(hexagons[SW][0])  # our S is south-west's NE (0)
                else:
                    vertices.append(corners[4])
                    hexagon.append(vertex_counter)
                    vertex_counter += 1
                # SE (5) corner could have been created by south-eastern neighbour
                if SE is not (False or -1):
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
            if self.alignment == 'square':
                raise Exception('PBCs not well '
                                'defined for a square arrangement')
                # if j == -1:
                #     i += int(self.ny / 2)
                #     j = self.ny - 1
                # if i == -1:
                #     i = int(self.ny / 2) - j
                #     j = self.ny - 1
                # if j == self.ny:
                #     j = 0
                #     i = int(self.ny / 2) - i
            if self.alignment == 'diagonal':
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

    def hexagon_corners(self, x, y, radius):
        """
        Returns the coordinates of the corners of the hexagonal with
        center at `x`, `y` and radius `radius`.

        """
        angle_deg = [60 * i + 30 for i in range(6)]
        angle_rad = [a * np.pi / 180 for a in angle_deg]
        return [(x + radius * np.cos(theta),
                 y + radius * np.sin(theta),
                 0) for theta in angle_rad]

    # -------------------------------------------------------------------------
    # Functions that return the number of neighbours at the n-th shell at
    # a lattice site with position indexes (i, j). For example,
    # the function _ngbs_first_shell returns the nearest neighbours

    def _ngbs_first_shell(self, i, j):

        # there can be periodicity along x and y axes, but not
        # along the third cubic axis (would be z), hence, we use
        # _index to disregard periodicity in the NW and SE directions.
        if self.alignment == 'diagonal':
            return [self.index(i + 1, j),       # right  east
                    self.index(i - 1, j),       # left   west
                    self.index(i, j + 1),       # up     north-east
                    self.index(i, j - 1),       # down   south-west
                    self._index(i - 1, j + 1),  #        north-west
                    self._index(i + 1, j - 1),  #        south-east
                    ]
        # For the square alignment, the neighbours position change
        # between odd and even rows
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 1, j),       # right  east
                        self.index(i - 1, j),       # left   west
                        self.index(i + 1, j + 1),   # up     north-east
                        self.index(i, j - 1),       # down   south-west
                        self.index(i, j + 1),       # .      north-west
                        self.index(i + 1, j - 1),   # .      south-east
                        ]
            else:
                return [self.index(i + 1, j),       # right  east
                        self.index(i - 1, j),       # left   west
                        self.index(i, j + 1),       # up     north-east
                        self.index(i - 1, j - 1),   # down   south-west
                        self.index(i - 1, j + 1),   # .      north-west
                        self.index(i, j - 1),       # .      south-east
                        ]

    def _ngbs_second_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 1, j + 1),       # north east
                    self.index(i - 1, j - 1),       # south west
                    self.index(i - 1, j + 2),       # up
                    self.index(i + 1, j - 2),       # down
                    self._index(i - 2, j + 1),      # north west
                    self._index(i + 2, j - 1),      # south east
                    ]
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 2, j + 1),   # north east
                        self.index(i - 1, j - 1),   # south west
                        self.index(i, j + 2),       # up
                        self.index(i, j - 2),       # down
                        self.index(i - 1, j + 1),   # north west
                        self.index(i + 2, j - 1),   # south east
                        ]
            else:
                return [self.index(i + 1, j + 1),   # north east
                        self.index(i - 2, j - 1),   # south west
                        self.index(i, j + 2),       # up
                        self.index(i, j - 2),       # down
                        self.index(i - 2, j + 1),   # north west
                        self.index(i + 1, j - 1)    # south east
                        ]

    def _ngbs_third_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 2, j),       # east
                    self.index(i - 2, j),       # west
                    self.index(i, j + 2),       # north east
                    self.index(i, j - 2),       # south west
                    self._index(i - 2, j + 2),  # north west
                    self._index(i + 2, j - 2)   # south east
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            return [self.index(i + 2, j),       # east
                    self.index(i - 2, j),       # west
                    self.index(i + 1, j + 2),   # north east
                    self.index(i - 1, j - 2),   # south west
                    self.index(i - 1, j + 2),   # north west
                    self.index(i + 1, j - 2)    # south east
                    ]

    def _ngbs_fourth_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 2, j + 1),   # north east 1
                    self.index(i - 2, j - 1),   # south west 1
                    self.index(i + 1, j + 2),   # north east 2
                    self.index(i - 1, j - 2),   # south west 2
                    self.index(i - 1, j + 3),   # north 1
                    self.index(i + 1, j - 3),   # south 1
                    self.index(i - 2, j + 3),   # north 2
                    self.index(i + 2, j - 3),   # south 2
                    self.index(i - 3, j + 2),   # north west 1
                    self.index(i + 3, j - 2),   # south east 1
                    self.index(i - 3, j + 1),   # north west 2
                    self.index(i + 3, j - 1)    # south east 2
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 3, j + 1),   # north east 1
                        self.index(i - 2, j - 1),   # south west 1
                        self.index(i + 2, j + 2),   # north east 2
                        self.index(i - 2, j - 2),   # south west 2
                        self.index(i + 1, j + 3),   # north 1
                        self.index(i, j - 3),       # south 1
                        self.index(i, j + 3),       # north 2
                        self.index(i + 1, j - 3),   # south 2
                        self.index(i - 2, j + 2),   # north west 1
                        self.index(i + 2, j - 2),   # south east 1
                        self.index(i - 2, j + 1),   # north west 2
                        self.index(i + 3, j - 1)    # south east 2
                        ]
            else:
                return [self.index(i + 2, j + 1),   # north east 1
                        self.index(i - 3, j - 1),   # south west 1
                        self.index(i + 2, j + 2),   # north east 2
                        self.index(i - 2, j - 2),   # south west 2
                        self.index(i, j + 3),       # north 1
                        self.index(i - 1, j - 3),   # south 1
                        self.index(i - 1, j + 3),   # north 2
                        self.index(i, j - 3),       # south 2
                        self.index(i - 2, j + 2),   # north west 1
                        self.index(i + 2, j - 2),   # south east 1
                        self.index(i - 3, j + 1),   # north west 2
                        self.index(i + 2, j - 1)    # south east 2
                        ]

    def _ngbs_fifth_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 3, j),           # east
                    self.index(i - 3, j),           # west
                    self.index(i, j + 3),           # north east
                    self.index(i, j - 3),           # south west
                    self._index(i - 3, j + 3),      # north west
                    self._index(i + 3, j - 3)       # south east
                    ]
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 3, j),       # east
                        self.index(i - 3, j),       # west
                        self.index(i + 2, j + 3),   # north-east
                        self.index(i - 1, j - 3),   # south-west
                        self.index(i - 1, j + 3),   # north-west
                        self.index(i + 2, j - 3)    # south-east
                        ]
            else:
                return [self.index(i + 3, j),       # east
                        self.index(i - 3, j),       # west
                        self.index(i + 1, j + 3),   # north-east
                        self.index(i - 2, j - 3),   # south-west
                        self.index(i - 2, j + 3),   # north-west
                        self.index(i + 1, j - 3)    # south-east
                        ]

    def _ngbs_sixth_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 2, j + 2),   # north east
                    self.index(i - 2, j - 2),   # south west
                    self.index(i - 2, j + 4),   # north
                    self.index(i + 2, j - 4),   # south
                    self._index(i - 4, j + 2),  # north west
                    self._index(i + 4, j - 2)   # south east
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            return [self.index(i + 3, j + 2),   # north east
                    self.index(i - 3, j - 2),   # south west
                    self.index(i, j + 4),       # north
                    self.index(i, j - 4),       # south
                    self.index(i - 3, j + 2),   # north west
                    self.index(i + 3, j - 2)    # south east
                    ]

    def _ngbs_seventh_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 3, j + 1),       # north east 1
                    self.index(i - 3, j - 1),       # south west 1
                    self.index(i + 1, j + 3),       # north east 2
                    self.index(i - 1, j - 3),       # south west 2
                    self.index(i - 1, j + 4),       # north 1
                    self.index(i + 1, j - 4),       # south 1
                    self.index(i - 3, j + 4),       # north 2
                    self.index(i + 3, j - 4),       # south 2
                    self.index(i - 4, j + 3),       # north west 1
                    self.index(i + 4, j - 3),       # south east 1
                    self.index(i - 4, j + 1),       # north west 2
                    self.index(i + 4, j - 1)        # south east 2
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 4, j + 1),       # north east 1
                        self.index(i - 3, j - 1),       # south west 1
                        self.index(i + 3, j + 3),       # north east 2
                        self.index(i - 2, j - 3),       # south west 2
                        self.index(i + 1, j + 4),       # north 1
                        self.index(i - 1, j - 4),       # south 1
                        self.index(i - 1, j + 4),       # north 2
                        self.index(i + 1, j - 4),       # south 2
                        self.index(i - 2, j + 3),       # north west 1
                        self.index(i + 3, j - 3),       # south east 1
                        self.index(i - 3, j + 1),       # north west 2
                        self.index(i + 4, j - 1)        # south east 2
                        ]
            else:
                return [self.index(i + 3, j + 1),       # north east 1
                        self.index(i - 4, j - 1),       # south west 1
                        self.index(i + 2, j + 3),       # north east 2
                        self.index(i - 3, j - 3),       # south west 2
                        self.index(i + 1, j + 4),       # north 1
                        self.index(i - 1, j - 4),       # south 1
                        self.index(i - 1, j + 4),       # north 2
                        self.index(i + 1, j - 4),       # south 2
                        self.index(i - 3, j + 3),       # north west 1
                        self.index(i + 2, j - 3),       # south east 1
                        self.index(i - 4, j + 1),       # north west 2
                        self.index(i + 3, j - 1)        # south east 2
                        ]

    def _ngbs_eigth_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 4, j),       # east
                    self.index(i - 4, j),       # west
                    self.index(i, j + 4),       # north east
                    self.index(i, j - 4),       # south west
                    self._index(i - 4, j + 4),  # north west
                    self._index(i + 4, j - 4)   # south east
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            return [self.index(i + 4, j),       # east
                    self.index(i - 4, j),       # west
                    self.index(i + 2, j + 4),   # north east
                    self.index(i - 2, j - 4),   # south west
                    self.index(i - 2, j + 4),   # north west
                    self.index(i + 2, j - 4)    # south east
                    ]

    def _ngbs_ninth_shell(self, i, j):

        if self.alignment == 'diagonal':
            return [self.index(i + 3, j + 2),       # north east 1
                    self.index(i - 3, j - 2),       # south west 1
                    self.index(i + 2, j + 3),       # north east 2
                    self.index(i - 2, j - 3),       # south west 2
                    self.index(i - 2, j + 5),       # north 1
                    self.index(i + 2, j - 5),       # south 1
                    self.index(i - 3, j + 5),       # north 2
                    self.index(i + 3, j - 5),       # south 2
                    self.index(i - 5, j + 3),       # north west 1
                    self.index(i + 5, j - 3),       # south east 1
                    self.index(i - 5, j + 2),       # north west 2
                    self.index(i + 5, j - 2)        # south east 2
                    ]
        # in this case, the index does not depend on the row
        elif self.alignment == 'square':
            if j % 2 == 0:
                return [self.index(i + 4, j + 2),       # north east 1
                        self.index(i - 4, j - 2),       # south west 1
                        self.index(i + 4, j + 3),       # north east 2
                        self.index(i - 3, j - 3),       # south west 2
                        self.index(i + 1, j + 5),       # north 1
                        self.index(i, j - 5),           # south 1
                        self.index(i, j + 5),           # north 2
                        self.index(i + 1, j - 5),       # south 2
                        self.index(i - 3, j + 3),       # north west 1
                        self.index(i + 4, j - 3),       # south east 1
                        self.index(i - 4, j + 2),       # north west 2
                        self.index(i + 4, j - 2)        # south east 2
                        ]
            else:
                return [self.index(i + 4, j + 2),       # north east 1
                        self.index(i - 4, j - 2),       # south west 1
                        self.index(i + 3, j + 3),       # north east 2
                        self.index(i - 4, j - 3),       # south west 2
                        self.index(i, j + 5),           # north 1
                        self.index(i - 1, j - 5),       # south 1
                        self.index(i - 1, j + 5),       # north 2
                        self.index(i, j - 5),           # south 2
                        self.index(i - 4, j + 3),       # north west 1
                        self.index(i + 3, j - 3),       # south east 1
                        self.index(i - 4, j + 2),       # north west 2
                        self.index(i + 4, j - 2)        # south east 2
                        ]
