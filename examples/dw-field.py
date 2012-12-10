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

def relax_system(mesh,mat):
    sim=Sim(mesh,T=0,mat=mat)
    sim.c=2
    
    exch=UniformExchange(mat.J,mu_s=mat.mu_s)
    sim.add(exch)

    anis=Anisotropy(mat.D,mu_s=mat.mu_s)
    sim.add(anis)
    
    sim.set_m(init_m)
    
    vs=VisualSpin(sim)
    vs.init()
    ts=np.linspace(0, 10, 2001)
    for t in ts:
        sim.run_until(t)
        vs.update()
        time.sleep(0.01)
        
    
    return sim.spin

def dw_motion(mesh,m0,mat,H0=1):
    sim=Sim(mesh,T=0,mat=mat)
    sim.c=3
    
    exch=UniformExchange(mat.J,mu_s=mat.mu_s)
    sim.add(exch)

    anis=Anisotropy(mat.D,mu_s=mat.mu_s)
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
    
    ni=Nickel()
    ni.alpha=0.5
    ni.a=1
    ni.b=1
    ni.c=1
    ni.unit_length=1
    ni.mu_s=1
    ni.D=1
    ni.J=0.2
    ni.gamma=10
    
    mesh=FDMesh(nx=30,ny=1,nz=1)
    mesh.set_material(ni)
    
    
    m0=relax_system(mesh,ni)
    print 'relax system done'
    #dw_motion(mesh,m0,ni)
    
    
    
