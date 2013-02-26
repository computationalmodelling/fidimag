from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Nickel
import numpy as np

def init_m(pos):
    x,y,z=pos
    if (x,y,z)==(1,2,3):
        return (1,2,3)
    elif z<1:
        return (0,0,-1)
    else:
        return (0,0,1)
    
def init_T(pos):
    return np.sum(pos)


def test_sim_init_m():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m((0,1,0))
    sim.spin.shape=(3,-1)
    spin_y=sim.spin[1]
    assert(spin_y.any()==1)
    


def test_sim_init_m_fun():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m(init_m,normalise=False)
    assert(sim.spin_at(1,2,3)==(1,2,3))


def test_sim_T_fun():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_T(init_T)
    assert(sim.T[0]==0)
    assert(sim.T[-1]==9)


if __name__=='__main__':
    #test_sim_init_m()
    #test_sim_init_m_fun()
    test_sim_T_fun()
