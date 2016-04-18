"""
example 2 in the nmag manual
from http://nmag.soton.ac.uk/nmag/0.2/manual/html/example2/doc.html

"""
import numpy as np
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Demag


mesh = CuboidMesh(3, 3, 10/3.0, 10, 10, 30, unit_length=1e-9)
sim = Sim(mesh, name="nmag2")
sim.Ms = 0.86e6
sim.alpha = 0.5
sim.set_m((1, 0, 1))
sim.add(UniformExchange(A=13e-12))
sim.add(Demag())

ts = np.linspace(0, 3e-10, 61)
for t in ts:
    sim.run_until(t)
