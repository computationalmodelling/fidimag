import numpy as np
#from pccp.pc import *
from pccp.pc import Sim
from pccp.pc import FDMesh
from pccp.pc import DMI
from pccp.pc import UniformExchange


def init_m(pos):
    x,y,z=pos
    if x<50:
        return (0,0,1)
    elif x>50-1:
        return (0,1,-1)
    else:
        return (0,1,0)


def relax_system(mesh):
        
    sim=Sim(mesh,name='relax',pbc='2d')
    sim.alpha=0.1

    J = 1
    exch = UniformExchange(J)
    sim.add(exch)

    dmi=DMI(0.1*J,direction=(1,0,0))
    sim.add(dmi)
    
    sim.set_m(init_m)
    
    ts=np.linspace(0, 5e-10, 101)
    for t in ts:
        print t,sim.spin_length()-1
        sim.run_until(t)
    
    sim.save_vtk()
    
    return sim.spin


                        
if __name__=='__main__':
    
    mesh=FDMesh(nx=20,ny=1,nz=1)
    
    m0=relax_system(mesh)
    print 'relax system done'
    #spin_wave(mesh,m0)
    
