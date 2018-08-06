import numpy as np
import fidimag
# PdFe on Ir(111) [PRL, 114(17):1-5, 2015]
Ms = 1.1e6
D  = -3.9e-3


XL = 15e-9
YL = XL
ZL = XL

nx =  2
ny = 1
nz = 1

def set_D(pos):
    x, y, z = pos
    if x < XL/2:
        return 0
    else:
        return -D

mesh = fidimag.common.CuboidMesh(nx=nx,ny=ny,nz=nz, dx=XL/nx, dy=YL/ny, dz=ZL/nz)

sim = fidimag.micro.Sim(mesh)
sim.set_Ms(Ms)


sim.add(fidimag.micro.DMI(set_D))
m0 = np.array([0, 0, 1])
sim.set_m(m0)

sim.compute_effective_field(0)
print(sim.driver.field.reshape(-1, 3))
