import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from pccp.pc import *
import time

def init_m(pos):
    x,y,z=pos
    if x<140*5:
        return (1,0,0)
    elif x>160*5:
        return (-1,0,0)
    else:
        return (0,1,0)
    

class Material(object):
    def __init__(self):
        self.a=5
        self.b=5
        self.c=5
        #Ms = 8e5, mu_s = Ms * v
        self.mu_s=1e-22
        
        #A =1e-11, J = 2a*A 
        self.J=1e-20
        self.Dx=0.005*self.J
        
        self.gamma=1.76e11
        self.alpha=0.01
        self.unit_length=1e-10
    

def relax_system(mesh):
    mat = Material()
    mesh.set_material(mat)
    
    sim=Sim(mesh,name='relax')
    sim.alpha=0.5
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)
    
    sim.set_m(init_m)
    
    ts=np.linspace(0, 1e-10, 101)
    for t in ts:
        sim.run_until(t)
        print t
        sim.save_vtk()
    
    np.save('m0.npy',sim.spin)
    

def save_plot():
    fig=plt.figure()
    data = np.load('m0.npy')
    data.shape=(3,-1)
    print data
    
    plt.plot(data[0],'-',label='Sx')
    plt.plot(data[1],'-',label='Sy')
    plt.plot(data[2],'-',label='Sz')

    plt.ylim([-1.2,1.2])
    plt.ylabel('S')
    plt.xlabel('x/a')
    plt.legend(loc=1)
    plt.grid()
    fig.savefig('mxyz.pdf')

if __name__=='__main__':
    
    mesh=FDMesh(nx=300)
    
    relax_system(mesh)
    print 'relax system done'
    save_plot()
    #spin_wave(mesh,m0)
    
