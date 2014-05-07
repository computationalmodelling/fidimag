import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange
from pccp.pc import Zeeman
from pccp.pc import Constant

const = Constant()

"""
If we only consider the exchange, dmi and external field, 
we don't have to consider the lattice constant a.
"""

def init_m(pos):
    x,y,z = pos
    
    x0,y0,r  = 70,25,10
    
    m1 = (0.05,0.01,-1)
    m2 = (0,0,1)
    
    if (x-x0)**2 + (y-y0)**2 < r**2:
        return m1
    else:
        return m2


def relax_system(mesh):
    
    sim=Sim(mesh,name='relax',pbc='2d')
    sim.set_options(gamma=const.gamma, k_B=const.k_B)
    sim.alpha = 0.5
    sim.mu_s = const.mu_s_1
    
    sim.set_m(init_m)
    #sim.set_m(random_m)
    #sim.set_m(np.load('m_10000.npy'))

    J = 1.0*const.k_B
    exch = UniformExchange(J)
    sim.add(exch)
    
    D = 0.5*J
    dmi = DMI(D)
    sim.add(dmi)
    
    Hz = 0.2*J/const.mu_s_1
    zeeman = Zeeman([0,0,Hz])
    sim.add(zeeman)
    
    ONE_DEGREE_PER_NS = 17453292.52
    
    sim.relax(dt=1e-13, stopping_dmdt=ONE_DEGREE_PER_NS, max_steps=1000, save_m_steps=100, save_vtk_steps=50)
    
    np.save('m0.npy',sim.spin)

def plot_mxyz(ts,mx,my,mz,me,name):
    fig=plt.figure()
    plt.plot(ts,mx,'^-',label='mx')
    plt.plot(ts,my,'.-',label='my')
    plt.plot(ts,mz,'o-',label='mz')
    plt.plot(ts,me,'>-',label='me')
    plt.legend()
    fig.savefig(name)


def temperature_test(T):
    ni = Nickel()
    ni.alpha=0.1
    ni.D = 1.35e-26
    ni.mu_s = 2.16e-23 
    ni.J = 6.16e-21
    
    (nx,ny,nz)=(24,24,24)

    mesh=FDMesh(nx=nx,ny=ny,nz=nz)
    

    sim=Sim(mesh,T=T, driver='sllg')
    sim.set_options(dt=1e-15,gamma=ni.gamma, k_B=ni.k_B)
    sim.mu_s = ni.mu_s

    exch=UniformExchange(ni.J)
    sim.add(exch)

    anis=Anisotropy(ni.D)
    sim.add(anis)

    #zeeman=Zeeman(1e2,(0,0,1))
    #sim.add(zeeman)

    #demag=Demag(mu_s=ni.mu_s)
    #sim.add(demag)

    sim.set_m((0,0.6,0.99))


    #vs=VisualSpin(sim)
    #vs.init()

    ts=np.linspace(0, 1e-11, 101)
    me=[]
    mx=[]
    my=[]
    mz=[]
    for t in ts:
        sim.run_until(t)
        #vs.update()
        av=sim.compute_average()
        mx.append(av[0])
        my.append(av[1])
        mz.append(av[2])
        me.append(np.sqrt(np.sum(av*av)))
        #time.sleep(0.001)
        print av,me[-1]
        sim.save_vtk()
    name='nx%d_ny%d_nz%d_T%g.png'%(nx,ny,nz,T)
    plot_mxyz(ts,mx,my,mz,me,name)


if __name__=='__main__':
    mesh=FDMesh(nx=150,ny=50,nz=1)
    relax_system(mesh)
