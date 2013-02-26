from pc import UniformExchange
from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Nickel
import numpy as np

def init_m(pos):
    x,y,z=pos
    return (x,y,z)
    

def test_exch_1d():
    mesh=FDMesh(nx=5,ny=1,nz=1)
    sim=Sim(mesh)
    exch=UniformExchange(1)
    sim.add(exch)
    
    sim.set_m(init_m,normalise=False)

    field=exch.compute_field()
    assert field[0]==1
    assert field[1]==2
    assert field[2]==4
    assert field[3]==6
    assert field[4]==3
    assert np.max(field[5:])==0

    
def test_exch_2d():
    mesh=FDMesh(nx=5,ny=2,nz=1)
    sim=Sim(mesh)
    exch=UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m,normalise=False)

    field=exch.compute_field()
    assert field[0]==1
    assert field[1]==2+1
    assert field[2]==1+2+3
    assert field[3]==2+3+4
    assert field[4]==3+4
    assert np.max(field[21:])==0


def test_exch_3d():
    mesh=FDMesh(nx=4,ny=3,nz=2)
    sim=Sim(mesh)
    exch=UniformExchange(1)
    sim.add(exch)

    sim.set_m(init_m,normalise=False)

    field=exch.compute_field()
    #print field
    assert field[0]==1
    assert field[1]==0+1+2+1
    assert field[2]==1+2+3+2
    assert field[3]==2+3+3
   
    assert field[4]==1
    assert field[5]==5
    assert field[6]==10
    assert field[7]==11
    


    
    
if __name__=='__main__':
    test_exch_1d()
    test_exch_2d()
    test_exch_3d()
