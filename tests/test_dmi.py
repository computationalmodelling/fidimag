from pc import DMI 
from pc import FDMesh
from pc import Sim
from pc import Nickel
import numpy as np


def test_dmi_1d():
    mesh=FDMesh(nx=2,ny=1,nz=1)
    
    sim=Sim(mesh)
    sim.set_m((1,0,0))
    
    dmi = DMI(D=1,direction=(0,0,1))
    sim.add(dmi)
    
    field=dmi.compute_field()

    expected=np.array([0,0,0,0,0,0])
    
    assert (field==expected).all()
   

    energy=dmi.compute_energy()
    assert energy==0


    
    
if __name__=='__main__':
    test_dmi_1d()
    
