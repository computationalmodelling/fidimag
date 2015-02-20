import numpy as np
from pccp.pc import *
import time

class Material(object):
    def __init__(self):
        self.a=5
        self.b=5
        self.c=5
        #Ms = 8e5, mu_s = Ms * v
        self.mu_s=1e-22
        
        #A =0.5e-11, J = 2a*A 
        self.J=1e-20
        self.Dx=0.003*self.J
        self.Dp=-0.02*self.J

        self.gamma=1.76e11
        self.alpha=0.01
        self.unit_length=1e-10

def init_m2(pos):
    return (0,0,1)
    x,y,z=pos
    if x<500:
        return (0,0,1)
    elif x>1500:
        return (0,0,-1)
    else:
        return (0,1,0)
    
def init_m(pos):
    return np.random.random(3)-0.5
    
def pin_fun(t,mesh,spin):
    n=mesh.nxyz
    
    spin[0]=0
    spin[n]=0
    spin[n+n]=1
    

def relax_system(mesh):
    
    mat = Material()
    mesh.set_material(mat)
    
    sim=Sim(mesh,name='dmi_2d')
    sim.alpha=0.1
    
    #sim.pin_fun=pin_fun
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    dmi=DMI(0.4*mat.J,direction=(1,1,0))
    sim.add(dmi)
    
    H0 = 2.5
    zeeman=Zeeman(H0,(0,0,1))
    sim.add(zeeman)
    
    #sim.set_m(init_m)
    sim.set_m(np.load('m0.npy'))
    
    ts=np.linspace(0, 1e-10, 101)
    for t in ts:
        print t
        sim.run_until(t)
        print sim.spin_length()
        sim.save_vtk()
    
    return sim.spin


                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=200,ny=30,nz=1)
    
    m0=relax_system(mesh)
    np.save('m0.npy',m0)
    print 'relax system done'
    #spin_wave(mesh,m0)
    
