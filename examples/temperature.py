import numpy as np
from pc import *
from pc.show_vector import VisualSpin
import time

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

def plot_mxyz(ts,mx,my,mz,me,name):
    fig=plt.figure()
    plt.plot(ts,mx,'^-',label='mx')
    plt.plot(ts,my,'.-',label='my')
    plt.plot(ts,mz,'o-',label='mz')
    plt.plot(ts,me,'>-',label='me')
    plt.legend()
    fig.savefig(name)


def temperature_test(T):
    ni=Nickel()
    ni.alpha=0.1
    ni.D*=0
    ni.mu_s*=1
    ni.J*=1
    (nx,ny,nz)=(10,10,10)

    mesh=FDMesh(nx=nx,ny=ny,nz=nz)
    mesh.set_material(ni)

    sim=Sim(mesh,T=T,mat=ni)

    exch=UniformExchange(ni.J)
    sim.add(exch)

    #anis=Anisotropy(ni.D,mu_s=ni.mu_s)
    #sim.add(anis)

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
    name='nx%d_ny%d_nz%d_T%g.png'%(nx,ny,nz,T)
    plot_mxyz(ts,mx,my,mz,me,name)






if __name__=='__main__':
    temperature_test(300)
