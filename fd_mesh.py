
class FDMesh():
    def __init__(self,dx=1,dy=1,dz=1,nx=10,ny=1,nz=1,unit_length=1):
        self.dx=dx
        self.dy=dy
        self.dz=dz
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nxyz=nx*ny*nz
        self.unit_length=unit_length
        


