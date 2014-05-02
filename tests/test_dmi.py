from pc import DMI 
from pc import FDMesh
from pc import Sim
from pc import Nickel
import numpy as np

def test_dmi_1d():
    mesh=FDMesh(nx=2,ny=1,nz=1)
    
    sim=Sim(mesh)
    sim.set_m((1,0,0))
    
    dmi = DMI(D=1)
    sim.add(dmi)
    
    field=dmi.compute_field()

    expected=np.array([0,0,0,0,0,0])
    
    assert (field==expected).all()
   
    energy=dmi.compute_energy()
    assert energy==0


def init_m(pos):
    x,y,z=pos
    if x<0.5:
        return (0,1,0)
    else:
        return (0,0,1)

def test_dmi_1d_field():
    
    mesh=FDMesh(nx=2,ny=1,nz=1)
    
    sim=Sim(mesh)
    sim.set_m(init_m)
    
    dmi = DMI(D=1.23)
    sim.add(dmi)
    
    field=dmi.compute_field()
    expected=np.array([0,0,-1,0,0,-1])*1.23
    
    assert (field==expected).all()
    
    energy=dmi.compute_energy()
    print energy
    assert energy==1.23
    
    
if __name__=='__main__':
    test_dmi_1d()
    test_dmi_1d_field()
    
