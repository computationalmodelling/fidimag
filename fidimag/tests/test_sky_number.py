from pc import Sim
from pc import FDMesh
from pc import DMI
from pc import UniformExchange
from pc import Zeeman


def init_m(pos):
    x = pos[0]
    y = pos[1]
    
    x0,y0,r  = 60,60,20
    
    m1 = (0.05,0.01,-1)
    m2 = (0,0,1)
    
    if (x-x0)**2 + (y-y0)**2 < r**2:
        return m1
    else:
        return m2


def test_skx_num():
    
    mesh = FDMesh(nx=120,ny=120,nz=1, pbc='xy')
    
    sim=Sim(mesh,name='skx_num')
    sim.set_tols(rtol=1e-6,atol=1e-6)
    sim.alpha = 1.0
    sim.gamma = 1.0
    sim.mu_s = 1.0
    
    sim.set_m(init_m)
    
    sim.do_procession = False

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)
    
    D = 0.09
    dmi = DMI(D)
    sim.add(dmi)
    
    zeeman = Zeeman([0,0,5e-3])
    sim.add(zeeman)
    
    sim.relax(dt=2.0, stopping_dmdt=1e-2, max_steps=1000, save_m_steps=None, save_vtk_steps=None)
    
    skn =  sim.skyrmion_number()
    print 'skx_number', skn
    assert skn > -1 and skn < -0.99 

                        
if __name__=='__main__':
    
    test_skx_num()
    
