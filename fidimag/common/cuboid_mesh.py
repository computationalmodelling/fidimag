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
from __future__ import print_function
from psutil import 	virtual_memory
import numpy as np
from textwrap import dedent
from six.moves import range


class CuboidMesh(object):

    def __init__(self, dx=1, dy=1, dz=1, nx=1, ny=1, nz=1, x0=0, y0=0, z0=0,
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
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.periodicity = periodicity

        self.Lx = dx * nx  # total size of mesh
        self.Ly = dy * ny
        self.Lz = dz * nz

        self.n = nx * ny * nz  # total number of cells
        self.nxy = nx * ny  # number of cells in the x-y plane
        self.size = (nx, ny, nz)
        self.check_size()

        self.mesh_type = "cuboid"
        self.unit_length = unit_length

        self.coordinates = self.init_coordinates()
        self.neighbours, self.next_neighbours = self.init_neighbours()
        self.grid = self.init_grid()  # for vtk export

    def __repr__(self):
        repres = dedent("""\
    Cuboid Mesh
    Dimensions = {} x {} x {}
    Discretisation = ({}, {}, {})
    (x0, y0, z0) = ({}, {}, {})
    Periodicity = {}
    Total number of cells = {}
    """)
        return repres.format(self.Lx, self.Ly, self.Lz,
                             self.dx, self.dy, self.dz,
                             self.x0, self.y0, self.z0,
                             self.periodicity,
                             self.n)

    def _repr_html_(self):
        repres = dedent("""\
    <h3>Cuboid Mesh:</h3>
    <b>Dimensions</b> = {} x {} x {} nm<br>
    <b>Discretisation</b> = ({}, {}, {})<br>
    <b>x0</b> = ({}, {}, {})<br>
    <b>Periodicity</b> = {}<br>
    <b>No. of Cells</b> = {}<br>
    """)
        return repres.format(self.dx, self.dy, self.dz,
                             self.nx, self.ny, self.nz,
                             self.x0, self.y0, self.z0,
                             self.periodicity,
                             self.n)

    def init_coordinates(self):
        coordinates = np.zeros((self.n, 3))
        for i in range(self.nz):
            for j in range(self.ny):
                for k in range(self.nx):
                    index = self.index(k, j, i)
                    r = (self.x0 + k * self.dx + self.dx / 2.0,
                         self.y0 + j * self.dy + self.dy / 2.0,
                         self.z0 + i * self.dz + self.dz / 2.0)
                    coordinates[index] = r
        return coordinates

    def init_grid(self, origin=(0, 0, 0)):
        """
        Compute the coordinates for each of the axes.

        """
        return (origin[0] + np.linspace(0, self.Lx, self.nx + 1),
                origin[1] + np.linspace(0, self.Ly, self.ny + 1),
                origin[2] + np.linspace(0, self.Lz, self.nz + 1))

    def init_neighbours(self):
        # array will have entry set to -1 for nonexisting neighbours
        # this way we get to use a 2d array which is convenient to use
        # in our C code instead of a list of lists
        connectivity = []
        connectivity_next = []
        for i in range(self.nz):
            for j in range(self.ny):
                for k in range(self.nx):
                    # Unique index for the current lattice site
                    cell = self._index(k, j, i)
                    neighbours = [
                        self.index(k - 1, j, i),  # left
                        self.index(k + 1, j, i),  # right
                        self.index(k, j - 1, i),  # behind
                        self.index(k, j + 1, i),  # in front
                        self.index(k, j, i - 1),  # under
                        self.index(k, j, i + 1),  # over
                    ]

                    # If one of the neighbours is the cell itself, we
                    # set its index to -1
                    # neighbours = [other if other != cell
                    #               else -1 for other in neighbours]

                    next_neighbours = [
                        self.index(k - 2, j, i),  # left
                        self.index(k + 2, j, i),  # right
                        self.index(k, j - 2, i),  # behind
                        self.index(k, j + 2, i),  # in front
                        self.index(k, j, i - 2),  # under
                        self.index(k, j, i + 2),  # over
                    ]

                    # July 1st, 2016 Weiwei: I think it's okay for a cell with
                    # its neighbour is itself if periodic boundary conditions
                    # are used. For example, if we only have one cell and
                    # enable periodic boundary condition in x-direction, then
                    # we got a rod.  therefore, I commented two lines below.
                    # no cell should be its own neighbour

                    connectivity.append(neighbours)
                    connectivity_next.append(next_neighbours)

        return (np.array(connectivity, dtype=np.int32),
                np.array(connectivity_next, dtype=np.int32)
                )

    def index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds. Handles periodic meshes.

        i, j, k are the positions in the x, y and z directions, respectively

        """
        if self.periodicity[0]:  # if mesh is periodic in x-direction
            if i < 0:            # then wrap the left side
                i += self.nx     # to the right
            elif i >= self.nx:   # and wrap the right side
                i -= self.nx     # to the left

        if self.periodicity[1]:
            if j < 0:
                j += self.ny
            elif j >= self.ny:
                j -= self.ny

        if self.periodicity[2]:
            if k < 0:
                k += self.nz
            elif k >= self.nz:
                k -= self.nz

        return self._index(i, j, k)

    def _index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or -1 if that cell would be out of bounds.

        """
        if i < 0 or j < 0 or k < 0 or k >= self.nz or j >= self.ny or i >= self.nx:
            return -1
        return k * self.nxy + j * self.nx + i

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

    def vertices(self, origin=(0, 0, 0)):
        """
        Compute the coordinates of all vertices.

        """
        grid = np.zeros(((self.nx + 1) * (self.ny + 1) * (self.nz + 1), 3))
        i = 0
        axes = self.axes(origin)
        for z in axes[2]:
            for y in axes[1]:
                for x in axes[0]:
                    grid[i] = (x, y, z)
                    i += 1
        return grid

    def check_size(self, system_memory_fake_for_testing=None):
        """
        Provide an estimate of how much system memory is needed and warn accordingly.

        """
        bytes_per_float_numpy = 8
        size_coordinates_bytes = self.nx * self.ny * \
            self.nz * 3 * bytes_per_float_numpy
        size_coordinates_GiB = size_coordinates_bytes / (1024. ** 3)

        if system_memory_fake_for_testing is None:
            mem = virtual_memory().total / (1024.0 ** 3)
        else:
            mem_GiB = system_memory_fake_for_testing

        if 2 * size_coordinates_GiB > mem_GiB:
            # print because no logging yet
            print("Warning! Size of mesh coordinates i {} GiB.".format(
                size_coordinates_GiB))
            print(
                "You have {} GiB system memory. Possible halt.".format(mem_GiB))
            return 1
        return 0
