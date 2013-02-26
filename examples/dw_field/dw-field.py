import numpy as np
from pc import *
import time


def init_m(pos):
    x,y,z=pos
    if x<8:
        return (1,0,0)
    elif x>12:
        return (-1,0,0)
    else:
        return (0,1,0)

def relax_system(mesh,mat):
    sim=Sim(mesh,T=0,mat=mat,name='relax')
    sim.c=2
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    anis=Anisotropy(mat.D)
    sim.add(anis)
    
    sim.set_m(init_m)
    
    ts=np.linspace(0, 100, 2001)
    sim.save_vtk()
    
    for t in ts:
        sim.run_until(t)
    
    sim.save_vtk()    
        
    return sim.spin

def dw_motion(mesh,m0,mat,H0=1):
    sim=Sim(mesh,T=0,mat=mat)
    sim.c=3
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    anis=Anisotropy(mat.D)
    sim.add(anis)
    
    zeeman=Zeeman(H0,(1,0,0))
    sim.add(zeeman)
    
    sim.set_m(m0)

    
    ts=np.linspace(0, 200, 2001)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
    
                        
if __name__=='__main__':
    
    ni=Nickel()
    ni.alpha=0.5
    ni.a=1
    ni.b=1
    ni.c=1
    ni.unit_length=1
    ni.mu_s=1
    ni.D=1
    ni.J=2
    ni.gamma=1
    
    mesh=FDMesh(nx=20,ny=5,nz=2)
    mesh.set_material(ni)
    
    
    m0=relax_system(mesh,ni)
    print 'relax system done'
    ni.alpha=0.05
    dw_motion(mesh,m0,ni)
    
    
    
