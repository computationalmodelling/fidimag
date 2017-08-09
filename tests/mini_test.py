import fidimag as f

mesh = f.common.CuboidMesh(nx=5, ny=1, nz=1)
sim = f.micro.Sim(mesh)

sim.set_m([0., 0., 1.])
