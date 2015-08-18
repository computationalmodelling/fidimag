import numpy as np

class FDMesh():
    """
    pbc should be a string contain 'x', 'y' or 'z'. '1d' and '2d' are also acceptable options.
    """
    def __init__(self, dx=1.0, dy=1.0, dz=1.0, nx=10, ny=1, nz=1, unit_length=1.0, pbc=None):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.nxy = nx*ny
        self.nyz = ny * nz
        self.nxyz = nx * ny * nz
        self.unit_length = unit_length
        
        self.compute_pos()

        self.pbc = pbc
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

                    tp = ((i + 0.5) * self.dx,
                          (j + 0.5) * self.dy,
                          (k + 0.5) * self.dz)

                    self.pos.append(tp)

    def _index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds.

        """
        if i < 0 or j < 0 or k < 0 or i >= self.nx or j >= self.ny or k >= self.nz:
            return -1
        return k * self.nxy + j * self.nx + i

    def index(self, i, j, k):
        """
        Returns the index for the cell with ordinals i, j, k
        or False if that cell would be out of bounds. Handles periodic meshes.
        """
        if self.xperiodic:  # if mesh is periodic in x-direction
            if i <= -1:          # then wrap the left side
                i = self.nx - 1  # to the right
            if i >= self.nx:     # and wrap the right side
                i = 0            # to the left
        
        if self.yperiodic:  
            if j <= -1:          
                j = self.ny - 1  
            if j >= self.ny:    
                j = 0

        if self.zperiodic:  
            if k <= -1:          
                k = self.nz - 1
            if k >= self.nz:    
                k = 0   

        return self._index(i, j, k)


    def init_neighbours(self):

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

    # only used for tests
    def pos_at(self, i, j, k):
        idx = self.index(i, j, k)
        return self.pos[idx]
