from pc import Anisotropy
from pc import FDMesh
from pc import Sim
from pc import Nickel

def init_m(pos):
    x,y,z=pos
    if (x,y,z)==(1,2,3):
        return (1,2,3)
    elif z<1:
        return (0,0,-1)
    else:
        return (0,0,1)




def test_sim_mesh():
    mesh=FDMesh(nx=3,ny=4,nz=5)
    sim=Sim(mesh)
    sim.set_m(init_m)
    sim.spin.shape=(3,-1)
    spin_z=sim.spin[2]
    spin_z.shape=(3,4,5)
    assert spin_z[1,2,3]==3




if __name__=='__main__':
    test_sim_mesh()
