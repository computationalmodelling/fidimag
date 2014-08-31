import numpy as np
from mesh import FDMesh
from sim import Sim
from save_vtk import SaveVTK
from materials import Nickel


def init_m(pos):
    x,y,z=pos
    return (0,np.cos(x*2e10),np.sin(x*2e10))

def test_save_vtk():
    T=1000
    ni=Nickel()
    ni.alpha=0.1

    mesh=FDMesh(nx=10,ny=2,nz=2)
    mesh.set_material(ni)
    print mesh.pos

    sim=Sim(mesh,T=T,mat=ni)
    sim.set_m(init_m)
    
    vtk = SaveVTK(mesh,sim.spin)
    vtk.save_vtk(sim.spin)
    


if __name__=='__main__':
    test_save_vtk()
