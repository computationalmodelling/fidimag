import numpy as np
from pc import *
import time

def init_T(pos):
    return 0.001*pos[0]

def init_m(pos):
    x,y,z=pos
    if x<300:
        return (1,0,0)
    elif x>350:
        return (-1,0,0)
    else:
        return (0,1,0)

def relax_system(mesh,mat):
    sim=Sim(mesh,T=0,mat=mat,name='relax')
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    anis=Anisotropy(mat.D)
    sim.add(anis)
    
    sim.set_m(init_m)
    
    ts=np.linspace(0, 5e-10, 21)
    sim.save_vtk()
    
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        
    return sim.spin

def dw_motion(mesh,m0,mat,H0=1):
    sim=Sim(mesh,T=0,mat=mat)
    sim.set_T(init_T)
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    anis=Anisotropy(mat.D)
    sim.add(anis)
    
    #zeeman=Zeeman(H0,(1,0,0))
    #sim.add(zeeman)
    
    sim.set_m(m0)

    
    ts=np.linspace(0, 1e-11, 201)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
    
                        
if __name__=='__main__':
    
    ni=Nickel()
    ni.J*=0.05
    
    mesh=FDMesh(nx=200,ny=1,nz=1)
    mesh.set_material(ni)
    
    m0=relax_system(mesh,ni)
    print 'relax system done'
    ni.alpha=0.05
    dw_motion(mesh,m0,ni)
    
    
    
