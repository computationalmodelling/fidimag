import numpy as np

from pc import FDMesh
from pc import Sim
from pc import Demag
from pc import UniformExchange

from util.oommf import compute_demag_field
from util.oommf import compute_exch_field


import clib

def test_oommf_coefficient():
    
    res =  clib.compute_Nxx(10,1,1,1,2,3)

    assert  abs(-0.000856757528962409449 - res) < 5e-15
    
    #print clib.compute_Nxx_asy(10,1,1,1,2,3)


def test_exch_field_oommf(A=1e-11, Ms=2.6e5):
    
    mesh = FDMesh(nx=10, ny=3, nz=2, dx=0.5, unit_length=1e-9)
    
    sim = Sim(mesh)
    sim.Ms = Ms
    
    exch = UniformExchange(A=A)
    sim.add(exch)
    
    def init_m(pos):
        
        x,y,z = pos
        
        return (np.sin(x)+y+2.3*z,np.cos(x)+y+1.3*z,0)
    
    sim.set_m(init_m)
    
    field = exch.compute_field()
    
    init_m0="""
    return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
    """
    field_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)
    
    print np.abs(max(field - field_oommf))
    assert np.abs(max(field - field_oommf))< 5e-7

def test_demag_field_oommf(Ms=6e5):
    mesh = FDMesh(nx=5,ny=2,nz=3, unit_length=1e-9)
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
    
    if { $x <=2e-9 } {
        return "1 0 0"
    } elseif { $x >= 4e-9 } {
        return "0 0 1"
    } else {
        return "0 1 0"
    }
    """
    
    field_oommf = compute_demag_field(mesh, Ms=Ms, init_m0=init_m0)
    assert np.max(abs(field - field_oommf)) < 2e-9
    
    
if __name__=='__main__':

    test_demag_field_oommf()
    
    test_exch_field_oommf()
    
