import os
import pyvtk
import numpy as np


class SaveVTK():
    def __init__(self,mesh,name='unnamed'):
        self.mesh=mesh
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.dx=mesh.dx
        self.dy=mesh.dy
        self.dz=mesh.dz
        self.name = '%s_vtks'%name
        xyz=np.array(mesh.pos)
        self.x=np.array(xyz[:,0],dtype='float32')
        self.y=np.array(xyz[:,1],dtype='float32')
        self.z=np.array(xyz[:,2],dtype='float32')
                
        #build a new index since we have used difference order 
        ids = [self.mesh.index(i,j,k) for k in range(self.nz) for j in range(self.ny) for i in range(self.nx)]
        self.ids = np.array(ids)
        
        self.pos = []
        for i in range(len(ids)):
            self.pos.append(self.mesh.pos[self.ids[i]])
    
    def save_vtk(self, m, step=0, vtkname='m'):
        
        if not os.path.exists(self.name):
            os.makedirs(self.name)
        
        pos=pyvtk.StructuredGrid([self.nx,self.ny,self.nz],self.pos)
        
        m.shape=(3,-1)
        data=pyvtk.PointData(pyvtk.Vectors(np.transpose(m)[self.ids],'m'))
        m.shape=(-1,)
        
        vtk = pyvtk.VtkData(pos,data,'spins')
                      
        vtk.tofile("%s/%s.%06d"%(self.name,vtkname,step),'binary')

    def save_vtk_scalar(self, skx_num, step=0, vtkname='skx'):
        
        folder_name = self.name + '_' + vtkname
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    
        pos=pyvtk.StructuredGrid([self.nx,self.ny,self.nz],self.pos)
    
        data=pyvtk.PointData(pyvtk.Scalars(skx_num[self.ids],'skx_num', lookup_table='default'))
    
        vtk = pyvtk.VtkData(pos,data,'skx_num')
    
        vtk.tofile("%s/%s.%06d"%(folder_name,vtkname,step),'binary')



