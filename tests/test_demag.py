from pc import FDMesh
from pc import Sim
from pc import Demag
from util.oommf import compute_demag_field
import numpy as np

import clib

def test_oommf_coefficient():
    
    res =  clib.compute_Nxx(10,1,1,1,2,3)

    assert  abs(-0.000856757528962409449 - res) < 5e-15
    
    #print clib.compute_Nxx_asy(10,1,1,1,2,3)

def test_demag_fft_exact():
    mesh = FDMesh(nx=5,ny=3,nz=4)
    sim = Sim(mesh)
    
    demag=Demag()
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
    #print fft,exact
    print np.max(np.abs(fft-exact))

    assert np.max(np.abs(fft-exact))<5e-22
    
def test_demag_fft_exact_oommf():
    mesh = FDMesh(nx=5,ny=3,nz=2)
    sim = Sim(mesh)
    
    demag=Demag(oommf=True)
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
    #print fft,exact
    print np.max(np.abs(fft-exact))

    assert np.max(np.abs(fft-exact))<1e-15


def test_demag_two_spin_xx():
    mesh = FDMesh(nx=2,ny=1,nz=1)
    sim = Sim(mesh)
    
    demag=Demag()
    sim.add(demag)
        
    sim.set_m((1,0,0))
    field=demag.compute_field()
    print field
    assert(field[0]==2e-7)
    assert(field[1]==2e-7)

def test_demag_field_oommf(Ms = 6e5 ):
    mesh = FDMesh(nx=5,ny=2,nz=3)
    sim = Sim(mesh)
    
    sim.Ms = Ms
    
    demag = Demag(oommf=True)
    sim.add(demag)
    
    def init_m(pos):
        
        x = pos[0]
        
        if x<=2:
            return (1,0,0)
        elif x>=4:
            return (0,0,1)
        else:
            return (0,1,0)
    
    sim.set_m(init_m)
    field = demag.compute_field()
    #exact=demag.compute_exact()
    
    init_m0="""
    
    if { $x <=2 } {
        return "1 0 0"
    } elseif { $x >= 4 } {
        return "0 0 1"
    } else {
        return "0 1 0"
    }
    """
    
    field_oommf = compute_demag_field(mesh, Ms=Ms, init_m0=init_m0)
    assert np.max(abs(field - field_oommf)) < 2e-9
    
    
if __name__=='__main__':
    #test_demag_fft_exact()
    #test_demag_two_spin_xx()
    
    #test_oommf_coefficient()
    #test_demag_fft_exact()
    #test_demag_fft_exact_oommf()
    test_demag_field_oommf()
    
