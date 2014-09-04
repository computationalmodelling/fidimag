

class FDMesh():
    def __init__(self,dx=1.0,dy=1.0,dz=1.0,nx=10,ny=1,nz=1,unit_length=1.0, pbc=None):
        self.dx=dx
        self.dy=dy
        self.dz=dz
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nyz=ny*nz
        self.nxyz=nx*ny*nz
        self.unit_length=unit_length
        self.compute_pos()
        
        self.pbc = pbc
        self.xperiodic = 0
        self.yperiodic = 0
        
        if pbc is None:
            pass
        elif pbc == 'x':
            self.xperiodic = 1
        elif pbc == 'y':
            self.yperiodic = 1
        elif pbc=='xy':
            self.xperiodic = 1
            self.yperiodic = 1
        else:
            raise Exception("options only can be None, 'x', 'y' or 'xy'.")
        
    def compute_pos(self):
        self.pos=[]
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):

                    tp=((i+0.5)*self.dx,
                        (j+0.5)*self.dy,
                        (k+0.5)*self.dz)
                    
                    self.pos.append(tp)
                    
    
    def index(self, i, j, k):
        idx = i*self.nyz + j*self.nz + k
        return idx

    #only used for tests
    def pos_at(self,i,j,k):
        idx = self.index(i,j,k)
        return self.pos[idx]
