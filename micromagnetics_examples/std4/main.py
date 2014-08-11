import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange, Demag
from pccp.pc import Zeeman, TimeZeeman
from pccp.pc import DataReader

mu0 = 4*np.pi*1e-7

def init_m(pos):
        
    x = pos[0]
        
    if x<=2:
        return (1,0,0)
    elif x>=4:
        return (0,0,1)
    else:
        return (0,1,0)

def relax_system(mesh):
    
    sim = Sim(mesh,name='relax')
    
    sim.set_options(rtol=1e-10,atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.0e5
    sim.do_procession = False
    
    #sim.set_m((1,0.25,0.1))
    sim.set_m(np.load('m0.npy'))
    
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)
    
    demag = Demag(oommf=True)
    sim.add(demag)
    
    mT = 795.7747154594767
            
    ONE_DEGREE_PER_NS = 17453292.52
    
    sim.relax(dt=1e-13, stopping_dmdt=0.0001*ONE_DEGREE_PER_NS, max_steps=5000, save_m_steps=100, save_vtk_steps=50)
    
    np.save('m0.npy',sim.spin)
    

def apply_field1(mesh):
    
    sim = Sim(mesh,name='dyn',pbc='2d')
    
    sim.set_options(rtol=1e-10,atol=1e-14)
    sim.alpha = 0.02
    sim.gamma = 2.211e5
    sim.Ms = 8.0e5
    
    sim.set_m(np.load('m0.npy'))
    
    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)
    
    demag = Demag(oommf=True)
    sim.add(demag)
    
    mT = 0.001/mu0
    print mT
    
    zeeman = Zeeman([-24.6*mT,4.3*mT,0],name='H')
    sim.add(zeeman, save_field=True)
        
    ts = np.linspace(0,1e-9,101)
    for t in ts:
        sim.run_until(t)
        print 'sim t=%g'%t
    
def deal_plot():
    data = DataReader('dyn.txt')
    ts = data['time']*1e9
    #mx = data['m_x']
    my = data['m_y']
    
    data2 = np.loadtxt('data.txt')
    
    ts2 = data2[:,0]
    my2 = data2[:,-2]
    
    plt.plot(ts, my, '.')
    plt.plot(ts2, my2, '--')
    
    #plt.legend()
    #plt.xlim([0, 0.012])
    #plt.ylim([-5, 100])
    plt.xlabel(r'Ts (ns)')
    #plt.ylabel('Susceptibility')
    plt.savefig('res.pdf')
    
                        
if __name__=='__main__':
    
    mesh = FDMesh(nx=200, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
    
    #relax_system(mesh)
    
    #apply_field1(mesh)
    deal_plot()
    
    
