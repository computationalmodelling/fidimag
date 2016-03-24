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

        self.mesh_type = "cuboid"
        self.unit_length = unit_length

        self.coordinates = self.init_coordinates()
        self.neighbours = self.init_neighbours()
        self.grid = self.init_grid()  # for vtk export

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
        for i in range(self.nz):
            for j in range(self.ny):
                for k in range(self.nx):
                    cell = self._index(k, j, i)
                    neighbours = [other for other in [
                        self.index(k - 1, j, i),  # left
                        self.index(k + 1, j, i),  # right
                        self.index(k, j - 1, i),  # behind
                        self.index(k, j + 1, i),  # in front
                        self.index(k, j, i - 1),  # under
                        self.index(k, j, i + 1),  # over
                    ]]
                    # no cell should be its own neighbour
                    neighbours = [other if other != cell
                                  else -1 for other in neighbours]
                    connectivity.append(neighbours)
        return np.array(connectivity, dtype=np.int32)

    def index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds. Handles periodic meshes.

        i, j, k are the positions in the x, y and z directions, respectively

        """
        if self.periodicity[0]:  # if mesh is periodic in x-direction
            if i == -1:          # then wrap the left side
                i = self.nx - 1  # to the right
            if i == self.nx:     # and wrap the right side
                i = 0            # to the left
        if self.periodicity[1]:
            if j == -1:
                j = self.ny - 1
            elif j == self.ny:
                j = 0
        if self.periodicity[2]:
            if k == -1:
                k = self.nz - 1
            if k == self.nz:
                k = 0
        return self._index(i, j, k)

    def _index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds.

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
