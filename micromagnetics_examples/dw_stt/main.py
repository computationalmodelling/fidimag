import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from micro import Sim
from micro import FDMesh
from micro import UniformExchange, Demag, DMI
from micro import Zeeman, TimeZeeman, UniaxialAnisotropy
from pc import DataReader

mu0 = 4*np.pi*1e-7


def init_m(pos):
        
    x = pos[0]
    
    if x<500:
        return (1,0,0)
    elif x<600:
        return (0,1,1)
    elif x<1200:
        return (-1,0,0)
    elif x<1300:
        return (0,1,1)
    else:
        return (1,0,0)
    
def init_m_1(pos):
        
    x = pos[0]
    
    if x<400:
        return (1,0,0)
    elif x<500:
        return (0,1,1)
    else:
        return (-1,0,0)


def relax_system(mesh):
    
    sim = Sim(mesh,name='relax')
    
    sim.set_tols(rtol=1e-8,atol=1e-10)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_procession = False
    
    sim.set_m(init_m_1)
    #sim.set_m(np.load('m0.npy'))
    
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    
    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)
    
    #dmi = DMI(D=8e-4)
    #sim.add(dmi)
    
    mT = 795.7747154594767
            
    ONE_DEGREE_PER_NS = 17453292.52
    
    sim.relax(dt=1e-14, stopping_dmdt=0.00001, max_steps=5000, save_m_steps=None, save_vtk_steps=None)
    
    np.save('m0.npy',sim.spin)
    
    
def deal_plot():

    data = np.load('m0.npy')
    
    data.shape = (3,-1)
    
    data = data[2]
    
    print data.shape
    
    data.shape = (80,80,2)
    
    mx = data[40,:,0]
    
    plt.plot(mx, '--', label='mz',dashes=(2,2))
    

    #plt.legend()
    #plt.xlim([0, 0.012])
    #plt.ylim([-5, 100])
    plt.xlabel(r'xs')
    #plt.ylabel('Susceptibility')
    plt.savefig('res.pdf')
    
def deal_plot2():

    
    all = []
    for i in range(500):
        data = np.load('dyn_npys/m_%d.npy'%i)
        data.shape = (3,-1)
        all.append(data[2][500])
    
    
    
    plt.plot(all, '--', label='mz')
    

    #plt.legend()
    #plt.xlim([0, 0.012])
    #plt.ylim([-5, 100])
    plt.xlabel(r't')
    #plt.ylabel('Susceptibility')
    plt.savefig('res.pdf')
    

def excite_system(mesh):
    
    sim = Sim(mesh,name='dyn', driver='llg_stt')
    
    sim.set_tols(rtol=1e-12,atol=1e-14)
    sim.alpha = 0.05
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    
    #sim.set_m(init_m)
    sim.set_m(np.load('m0.npy'))
    
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)
    
    anis = UniaxialAnisotropy(5e4)
    sim.add(anis)
    
    #dmi = DMI(D=8e-4)
    #sim.add(dmi)
    
    sim.jx = -1e12
    sim.beta = 1
    
    ts = np.linspace(0, 5e-9, 501)
    
    for t in ts:
        print 'time', t
        sim.run_until(t)
        sim.save_vtk()
        sim.save_m()
                        
if __name__=='__main__':
    
    mesh = FDMesh(nx=1000, ny=1, nz=1, dx=2, dy=2, dz=2.0, unit_length=1e-9)
    
    deal_plot2()
     
    #relax_system(mesh)
    
    #deal_plot()
    
    #excite_system(mesh)
    
    #apply_field1(mesh)
    #deal_plot()
    
    
