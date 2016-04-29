from fidimag.atomistic import Sim
from fidimag.common.cuboid_mesh import CuboidMesh
from fidimag.atomistic import UniformExchange, Zeeman
import fidimag.common.constant as const

mesh = CuboidMesh(nx=1, ny=1, dx=1, dy=1)

sim = Sim(mesh, name='relax_sk')
sim.gamma = const.gamma
sim.set_m((1, 0, 0))
sim.add(Zeeman((0, 0, 25.)))

sim.run_until(1e-11)
sim.set_tols(rtol=1e-10, atol=1e-12)
sim.run_until(2e-11)
