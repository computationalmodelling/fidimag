import numpy as np
from pyvtk import *
import pyvtk

class SaveVTK():
    def __init__(self,mesh,m,name='unnamed'):
        self.mesh=mesh
        self.m=m
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.dx=mesh.dx
        self.dy=mesh.dy
        self.dz=mesh.dz
        self.name=name
        self.index=0
        xyz=np.array(mesh.pos)
        self.x=np.array(xyz[:,0],dtype='float32')
        self.y=np.array(xyz[:,1],dtype='float32')
        self.z=np.array(xyz[:,2],dtype='float32')
        

        
    def save_vtk(self,m):
        
        pos=pyvtk.StructuredGrid([self.nx,self.ny,self.nz],self.mesh.pos)
        
        m.shape=(3,-1)
        data=pyvtk.PointData(pyvtk.Vectors(np.transpose(m),'m'))
        m.shape=(-1,)
        
        vtk = pyvtk.VtkData(pos,data,'spins')
                      
        vtk.tofile("%s.%06d"%(self.name,self.index))
        
        self.index+=1