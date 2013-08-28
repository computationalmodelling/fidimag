from materials import UnitMaterial

class FDMesh():
    def __init__(self,dx=1.0,dy=1.0,dz=1.0,nx=10,ny=1,nz=1,unit_length=1.0):
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
        
        self.mat = UnitMaterial()
        self.mat.a = dx
        self.mat.b = dy
        self.mat.c = dz
        self.mat.unit_length = unit_length
        
    def set_material(self,mat):
        self.unit_length=mat.unit_length
        self.dx=mat.a
        self.dy=mat.b
        self.dz=mat.c
        self.compute_pos()
        self.mat = mat
        
    def compute_pos(self):
        self.pos=[]
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):

                    tp=(i*self.dx,
                        j*self.dy,
                        k*self.dz)
                    
                    self.pos.append(tp)
                    
    
    def index(self, i, j, k):
        idx = i*self.nyz + j*self.nz + k
        return idx

    #only used for tests
    def pos_at(self,i,j,k):
        idx=self.index(i,j,k)
        return self.pos[idx]
