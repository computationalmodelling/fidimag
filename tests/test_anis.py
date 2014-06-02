from pc import Anisotropy
from pc import FDMesh
from pc import Sim
import numpy as np


def test_anis():
    mesh=FDMesh(nx=5,ny=3,nz=2)
    spin=np.zeros(90)
    anis=Anisotropy([1,0,0])
    anis.setup(mesh,spin,np.array([1]))
    field=anis.compute_field()
    assert len(mesh.pos)==5*3*2
    assert np.max(field)==0
    spin[0]=99
    field=anis.compute_field()
    assert field[0]==2*99


    
    
if __name__=='__main__':
    test_anis()
    
