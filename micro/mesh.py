import numpy as np

class FDMesh():
    def __init__(self,dx=1.0,dy=1.0,dz=1.0,nx=10,ny=1,nz=1,unit_length=1.0, x0=0, y0=0, z0=0, pbc=None):
        """
        pbc could be None, 'x', 'y' or 'xy'.
        """
        self.dx=dx
        self.dy=dy
        self.dz=dz
        self.dx_real = dx*unit_length
        self.dy_real = dy*unit_length
        self.dz_real = dz*unit_length
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nyz=ny*nz
        self.nxyz=nx*ny*nz
        self.unit_length=unit_length

        self.pbc = pbc
        
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        
        self.compute_pos()
        
        self.cellsize=dx*dy*dz*unit_length**3
        
        self.xperiodic = 0
        self.yperiodic = 0
    
        if pbc == 'x':
            self.xperiodic = 1
        elif pbc == 'y':
            self.yperiodic = 1
        elif pbc=='xy':
            self.xperiodic = 1
            self.yperiodic = 1
        else:
            raise Exception("Only options 'x', 'y' or 'xy' are acceptable!")
        
    def compute_pos(self):
        self.pos=[]
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):

                    tp=((i+0.5)*self.dx+self.x0,
                        (j+0.5)*self.dy+self.y0,
                        (k+0.5)*self.dz+self.z0)
                    
                    self.pos.append(tp)
                    
    
    def index(self, i, j, k):
        idx = i*self.nyz + j*self.nz + k
        return idx
    
    def index_z(self, k=0):
        ids = [self.index(i,j,k) for i in range(self.nx) for j in range(self.ny)]
        return np.array(ids)

    #only used for tests
    def pos_at(self,i,j,k):
        idx = self.index(i,j,k)
        return self.pos[idx]
