import numpy as np
from pccp import *
import time


def init_m(pos):
    x,y,z=pos
    if x<-2:
        return (1,0,0)
    elif x>2:
        return (-1,0,0)
    else:
        return (0,1,0)

def relax_system(mesh):
    sim=Sim(mesh)
    sim.alpha=0.5
    sim.gamma=5
    
    exch=UniformExchange(1)
    sim.add(exch)

    anis=Anisotropy(0.4)
    sim.add(anis)
    
    sim.set_m(init_m)
    
    ts=np.linspace(0, 20, 2001)
    for t in ts:
        sim.run_until(t)
        
    
    return sim.spin

def dw_motion(mesh,m0,H0=1):
    sim=Sim(mesh)
    sim.alpha=0.1
    sim.gamma=5
    
    exch=UniformExchange(1)
    sim.add(exch)

    anis=Anisotropy(0.4)
    sim.add(anis)
    
    zeeman=Zeeman(H0,(1,0,0))
    sim.add(zeeman)
    
    sim.set_m(m0)
    vs=VisualSpin(sim)
    vs.init()

    
    ts=np.linspace(0, 20, 2001)
    for t in ts:
        sim.run_until(t)
        vs.update()
        time.sleep(0.01)
    
                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=30,ny=5,nz=2)
    
    m0=relax_system(mesh)
    print 'relax system done'
    dw_motion(mesh,m0)
    
    
    
