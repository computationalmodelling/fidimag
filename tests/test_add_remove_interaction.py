from fidimag.common import CuboidMesh, constant
from fidimag.micro import UniformExchange, Sim, Zeeman
import numpy as np

def test_add_remove_interaction_simple():
    mesh = CuboidMesh(nx=10, ny=10, nz=10, unit_length=1e-9)
    name = 'test_add_remove_intn_simple'
    sim = Sim(mesh, name=name)
    sim.set_m(lambda pos: (0, 0, 1))
    sim.set_Ms(5.8e5)
    exch = UniformExchange(A=1e-11, name='Exchange')
    zee = Zeeman((0, 0, 0.05 / constant.mu_0), name='Zeeman')
    sim.add(exch)
    sim.add(zee)
    sim.driver.run_until(1e-9)
    sim.remove('Zeeman')
    sim.driver.run_until(2e-9)
    f = open(name + '.txt')
    lines = f.read().split('\n')
    headers = lines[0].split()
    first_data = lines[2].split()
    last_data = lines[2].split()
    # Find the position in the data table
    position = headers.index('E_Zeeman')
    assert np.abs(float(last_data[position])) < 1e-15

if __name__ == '__main__':
    test_add_remove_interaction_simple()