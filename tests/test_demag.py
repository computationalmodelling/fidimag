from pc import FDMesh
from pc import Sim
from pc import Demag
import numpy as np

def test_demag_fft_exact():
    mesh = FDMesh(nx=5,ny=3,nz=2)
    sim = Sim(mesh)
    
    demag=Demag(mu_s=1e3)
    sim.add(demag)
    
    def init_m(pos):
        x,y,z=pos
        if x<=2:
            return (1,0,0)
        elif x>=4:
            return (0,0,1)
        else:
            return (0,1,0)
    
    sim.set_m(init_m)
    fft=demag.compute_field()
    exact=demag.compute_exact()
    #print np.max(np.abs(fft-exact))

    assert np.max(np.abs(fft-exact))<5e-16


def test_demag_two_spin_xx():
    mesh = FDMesh(nx=2,ny=1,nz=1)
    sim = Sim(mesh)
    
    demag=Demag(mu_s=1e3)
    sim.add(demag)
        
    sim.set_m((1,0,0))
    field=demag.compute_field()
    print field
    assert(field[0]==0.2)
    assert(field[1]==0.2)


    
    
if __name__=='__main__':
    test_demag_fft_exact()
    test_demag_two_spin_xx()
    
