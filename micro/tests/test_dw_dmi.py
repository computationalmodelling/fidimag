import numpy as np

from micro import FDMesh
from micro import Sim
from micro import UniformExchange
from micro import UniaxialAnisotropy
from micro import DMI

import matplotlib.pyplot as plt

mesh = FDMesh(dx=2, nx=150, x0=-150, unit_length=1e-9)

def m_init_dw(pos):
    
    x = pos[0]
    
    if x < -10:
        return (1,0,0)
    elif  x > 10:
        return (-1,0,0)
    else:
        return (0,1,0)

def analytical(xs, A = 1.3e-11, D = 4e-4, K = 8e4):
        
    delta = np.sqrt(A/(K-D*D/(4*A)))*1e9
    
    phi = D/(2*A)*xs*1e-9
    
    mx = - np.tanh(xs/delta)
    my = 1.0/np.cosh(xs/delta)*np.cos(phi)
    mz = 1.0/np.cosh(xs/delta)*np.sin(phi)
    return mx,my,mz
    
def save_plot():
    fig=plt.figure()
    mxyz=np.load('relaxed.npy')
    mxyz.shape=(3,-1)
    xs = np.array([p[0] for p in mesh.pos])
    plt.plot(xs,mxyz[0],'.',label='mx')
    plt.plot(xs,mxyz[1],'.',label='my')
    plt.plot(xs,mxyz[2],'.',label='mz')
    
    mx,my,mz = analytical(xs)
    plt.plot(xs,mx,'-')
    plt.plot(xs,my,'-')
    plt.plot(xs,mz,'-')
    
    plt.ylabel('mxyz')
    plt.xlabel('x (nm)')
    plt.legend(loc=1)
    plt.grid()
    plt.xlim([-150,150])
    plt.ylim([-1.2,1.2])
    fig.savefig('dw1.pdf')

def relax_system(mesh=mesh):
    
    Ms = 8.0e5
    sim = Sim(mesh, name = 'relax')
    
    sim.set_options(rtol=1e-8,atol=1e-12)
    sim.Ms = Ms
    sim.alpha = 0.5
    #sim.do_procession = False
    
    sim.set_m(m_init_dw)
    
    A = 1.3e-11
    D = 4e-4
    Kx = 8e4
    Kp = -6e5
    
    sim.add(UniformExchange(A))
    sim.add(DMI(D))
    sim.add(UniaxialAnisotropy(Kx,axis=[1,0,0], name='Kx'))
    
    
    sim.relax(stopping_dmdt=0.01)
    #df.plot(sim.llg._m)
    np.save('relaxed.npy',sim.spin)
    #sim.save_vtk()
    save_plot()
    #df.interactive()
    
    
    
if __name__ == '__main__':
    
    relax_system()
    save_plot()
    #excite_system()
    #excite_system()

