from __future__ import division

import numpy as np
from visual import *
#import visual as vis

class VisualSpin(object):
    def __init__(self,sim):
        self.sim=sim
        mesh=sim.mesh
        self.mesh=mesh
        self.spin=sim.spin
        self.dx=mesh.dx
        self.dy=mesh.dy
        self.dz=mesh.dz
        self.nx=mesh.nx
        self.ny=mesh.ny
        self.nz=mesh.nz
        self.nxyz=mesh.nxyz
        
    
    def init(self):
        self.vspins=[]
        tmp=np.reshape(self.spin, (self.nxyz,3), order='F')
        nxy=self.nx*self.ny
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    pos=(i*self.dx,j*self.dy,k*self.dz)
                    id=k*nxy+j*self.ny+i
                    a=arrow(pos=pos, axis=vector(tmp[id]), color=color.yellow)
                    self.vspins.append(a)
        
    def update(self):
        tmp=np.reshape(self.spin, (self.nxyz,3), order='F')
        for i in range(len(tmp)):
            self.vspins[i].axis=vector(tmp[i])

