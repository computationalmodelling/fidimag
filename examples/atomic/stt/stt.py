import numpy as np
from pc import *

class Material(object):
    def __init__(self):
        self.a=5
        self.b=5
        self.c=5
        #Ms = 8e5, mu_s = Ms * v
        self.mu_s=1e-22
        
        #A =0.5e-11, J = 2a*A
        self.J=1e-20
        self.Dx=0.005*self.J
        self.Dp=-0.02*self.J
        
        self.gamma=1.76e11
        self.alpha=0.01
        self.unit_length=1e-10

def init_m(pos):
    x,y,z=pos
    if x<250*5:
        return (1,0,0)
    elif x>270*5:
        return (-1,0,0)
    else:
        return (0,1,0)


def relax_system(mesh):
    mat = Material()
    mesh.set_material(mat)
    sim=Sim(mesh)
    
    sim.set_m(init_m)
    exch=UniformExchange(mat.J)
    sim.add(exch)
    
    Kx=Anisotropy(mat.Dx, direction=(1,0,0),name='Dx')
    sim.add(Kx)
    
    Kp=Anisotropy(mat.Dp, direction=(0,0,1),name='Dp')
    sim.add(Kp)
    
    sim.alpha = 0.5
    
    ts = np.linspace(0, 2e-10, 101)
    for t in ts:
        sim.run_until(t)
    sim.save_vtk()
    
    np.save('m0.npy',sim.spin)

def dynamic(mesh):
    mat = Material()
    mesh.set_material(mat)
    
    sim=Sim(mesh, driver='llg_stt',name='stt')
    
    m0 = np.load('m0.npy')
    sim.set_m(m0)
    exch=UniformExchange(mat.J)
    sim.add(exch)
    
    Kx=Anisotropy(mat.Dx, direction=(1,0,0),name='Dx')
    sim.add(Kx)
    
    Kp=Anisotropy(mat.Dp, direction=(0,0,1),name='Dp')
    sim.add(Kp)
    
    sim.jx = 2e12
    sim.alpha=0.01
    sim.beta=0.02
    
    ts = np.linspace(0, 1e-9, 101)
    for t in ts:
        sim.run_until(t)
        sim.save_vtk()
        print t



if __name__=='__main__':
    mesh=FDMesh(nx=500,ny=1,nz=1)
    relax_system(mesh)
    dynamic(mesh)
