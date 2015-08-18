

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

    def index(self, i, j, k):

    	idx = k * self.nxy + j * self.nx + i
        
        return idx

    def init_neighbours(self):
        neighbours_x = []
        neighbours_y = []
        neighbours_z = []
        neighbours = []
        for i in xrange(self.nx):
            for j in xrange(self.ny):
                for k in xrange(self.nz):
                    index = k * self.nx*self.ny + j * self.nx + i

                    ngx = []
                    if i > 0 or self.xperiodic:
                        idx = index - self.nyz
                        if i == 0:
                            idx += self.nxyz
                        ngx.append(idx)

                    if i < self.nx - 1 or self.xperiodic:
                        idx = index + self.nyz
                        if i == self.nx-1:
                            idx -= self.nxyz
                        ngx.append(idx)

                    neighbours_x.append(ngx)

                    ngy = []
                    if j > 0 or self.yperiodic:
                        idy = index - self.nz
                        if j == 0: 
                            idy += self.nyz
                        ngy.append(idy)

                    if j < self.ny - 1 or self.yperiodic:
                        idy = index + self.nz
                        if j == self.ny-1:
                            idy -= self.nyz
                        ngy.append(idy)
                    neighbours_y.append(ngy)

                    ngz = []
                    if k>0:
                        ngz.append(k-1)
                    if k< self.nz-1:
                        ngz.append(k+1)

                    neighbours_z.append(ngz)
                    
                    neighbours.append(ngx+ngy+ngz)

        self.neighbours_x = neighbours_x
        self.neighbours_y = neighbours_y
        self.neighbours_z = neighbours_z
        self.neighbours = neighbours

    # only used for tests
    def pos_at(self, i, j, k):
        idx = self.index(i, j, k)
        return self.pos[idx]
