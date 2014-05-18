import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from pc import UniformExchange
from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Zeeman
from pc import UnitMaterial
from pc import DataReader
import numpy as np

class UnitMaterial(object):
    def __init__(self):
        self.a=1
        self.b=1
        self.c=1
        self.J=1.0
        self.Dx=0.5
        self.mu_s=1.0
        self.gamma=2*np.pi
        self.alpha=0.01
        self.unit_length=1.0

def init_m(pos):
    x,y,z=pos
    if x<12:
        return (1,0,0)
    elif x>18:
        return (-1,0,0)
    else:
        return (0,1,0)

def pin_fun(pos):
    if pos[0]==0:
        return 1
    else:
        return 0    


def relax_system(mesh,Dx=0.005,Dp=0.01):
    
    mat = UnitMaterial()
    
    sim = Sim(mesh,name='test_energy')
    sim.set_options(rtol=1e-10,atol=1e-12)
    
    sim.alpha = mat.alpha
    sim.gamma = mat.gamma
    sim.pins = pin_fun
    
    
    exch=UniformExchange(mat.J)
    sim.add(exch)

    anis=Anisotropy(Dx, direction=(1,0,0), name='Dx')
    sim.add(anis)
    
    anis2=Anisotropy(-Dp, direction=(0,0,1), name='Dp')
    sim.add(anis2)
    
    sim.set_m((1,1,1))
    
    T = 100
    ts=np.linspace(0, T, 201)
    for t in ts:
        sim.save_vtk()
        sim.run_until(t)
    

    sim.save_vtk()
    np.save('m0.npy', sim.spin)


def save_plot():
    fig=plt.figure()
        
    data = DataReader('relax.txt')
    ts = data['time']
    
    #plt.plot(ts, data['E_Dp'], label='E_Dp')
    #plt.plot(ts, data['E_Dp'], label='E_Dp')
    #plt.plot(ts, data['E_Dx'], label='E_Dx')
    plt.plot(ts, data['E_exch'], label='E_exch')
    plt.plot(ts, data['E_total'], label='E_total')
    
    
    plt.ylabel('Energy (J)')
    plt.xlabel('Time (s)')
    plt.legend(loc=1)
    plt.grid()
    #plt.ylim([-1.2,1.2])
    fig.savefig('energy.pdf')
    
def test_energy(do_plot=False):
    mesh=FDMesh(nx=30)
    relax_system(mesh)
    
    if do_plot:
        save_plot()
    
    data = DataReader('relax.txt')
    energy = data['E_total']
    
    for i in range(len(energy)-1):
        assert energy[i]>energy[i+1]

    
                        
if __name__=='__main__':
    
    test_energy(do_plot=True)
