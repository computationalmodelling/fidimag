"""
example 2 in the nmag manual
from http://nmag.soton.ac.uk/nmag/0.2/manual/html/example2/doc.html

"""
import numpy as np
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Demag

mesh = CuboidMesh(3, 3, 10/3.0, 10, 10, 30, unit_length=1e-9)


def run(integrator, jacobian):
    name = "sim_" + integrator
    if integrator == "sundials":
        name += "_J1" if jacobian else "_J0"
    sim = Sim(mesh, name, integrator=integrator, driver='llg', use_jac=jacobian)
    sim.Ms = 0.86e6
    sim.driver.alpha = 0.5
    sim.set_m((1, 0, 1))
    sim.add(UniformExchange(A=13e-12))
    sim.add(Demag())

    ts = np.linspace(0, 3e-10, 61)
    for t in ts:
        print(f'Running until t = {t}')
        sim.driver.run_until(t)

if __name__ == "__main__":
    run("sundials", False)
    run("sundials_diag", False)
    run("sundials", True)
