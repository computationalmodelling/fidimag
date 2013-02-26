class FDMesh():
    def __init__(self,dx=1.0,dy=1.0,dz=1.0,nx=10,ny=1,nz=1,unit_length=1.0):
        self.dx=dx
        self.dy=dy
        self.dz=dz
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nxy=nx*ny
        self.nxyz=nx*ny*nz
        self.unit_length=unit_length
        self.compute_pos()
        
    def set_material(self,mat):
        self.unit_length=mat.unit_length
        self.dx=mat.a/self.unit_length
        self.dy=mat.b/self.unit_length
        self.dz=mat.c/self.unit_length
        self.compute_pos()
        
    def compute_pos(self):
        self.pos=[]
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                
                    tp=(i*self.dx,
                        j*self.dy,
                        k*self.dz)
                    #index=k*self.nxy+j*self.ny+i
                    self.pos.append(tp)

    #only used for tests
    def pos_at(self,i,j,k):
        index=k*self.nxy+j*self.nx+i
        return self.pos[index]