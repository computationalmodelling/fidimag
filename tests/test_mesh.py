from pc import FDMesh
from pc import Nickel

def test_mesh1():
    mesh=FDMesh(nx=5,ny=3,nz=2,dx=0.23,dy=0.41)
    assert len(mesh.pos)==5*3*2
    assert mesh.nyz==6
    assert mesh.pos_at(0,0,0)==(0,0,0)
    assert mesh.pos_at(3,2,1)==(3*0.23,2*0.41,1)

if __name__=='__main__':
    test_mesh1()
