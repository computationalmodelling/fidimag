import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange
from pccp.pc import Demag
from pccp.pc import Anisotropy
from pccp.pc import Constant
from pccp.pc.save_vtk import SaveVTK

from pccp.pc.batch_task import BatchTasks

const = Constant()


def mu_s(pos):
    x,y,z = pos
    x0,y0,r  = 60*0.5,60*0.5,30*0.5
    
    if (x-x0)**2 + (y-y0)**2 <= r**2:
        return const.mu_s_1
    else:
        return 0

def init_m(pos):
    x,y,z = pos
    
    return (0,0,1)
    

def random_m(pos):
    return np.random.random(3)-0.5

def plot_f(mesh,field, mu_s_inv):
    fig = plt.figure(figsize=(6,4))
    
    field.shape = (3,-1)
    
    xs = range(121)
    fz = []
    for i in xs:
        j=mesh.index(i,60,0)
        if mu_s_inv[j]>0:
            fz.append(0)
        else:
            fz.append(field[2][j])

    
    plt.plot(xs,fz)
    plt.savefig('f.pdf')
    
def relax_system():
    
    mesh = FDMesh(nx=121,ny=121,nz=1,dx=0.5,dy=0.5, unit_length=1e-9)
    
    sim=Sim(mesh,name='relax_skx')
    sim.set_options(gamma=const.gamma, k_B=const.k_B)
    
    sim.alpha = 1.0

    sim.mu_s = mu_s
    
    sim.set_m(init_m)
    
    demag = Demag()
    
    sim.add(demag)
    
    field = demag.compute_field()
    
    print field
    
    vtk = SaveVTK(mesh, field, name='demag')
    vtk.save_vtk(field)
    
    sim.save_vtk()
    
    plot_f(mesh,field, sim.mu_s_inv)
    
    
    #np.save('m0.npy',sim.spin)



                        
if __name__=='__main__':
    
    np.random.seed(3)
    
    relax_system()
    
