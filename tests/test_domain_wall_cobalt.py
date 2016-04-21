import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, UniaxialAnisotropy

# material paramters
Ms_Co = 1400e3  # A/m
K1_Co = 520e3  # A/m
A_Co = 30e-12  # J/m

LENGTH = 100
UNIT_LENGTH = 1e-9
NODE_COUNT = 100


def initial_m(r):
    x = r[0]
    mz = 1.0 - 2.0 * x / LENGTH
    my = math.sqrt(1 - mz * mz)
    return [0, my, mz]


def reference_mz(x):
    # analytical solution for the relaxed mz
    return math.cos(math.pi / 2 + math.atan(math.sinh(((x - LENGTH / 2) * UNIT_LENGTH) / math.sqrt(A_Co / K1_Co))))


def setup_domain_wall_cobalt(node_count=NODE_COUNT, A=A_Co, Ms=Ms_Co, K1=K1_Co, length=LENGTH, do_precession=True, unit_length=UNIT_LENGTH):
    a = length / node_count  # cell size
    mesh = CuboidMesh(dx=a, dy=a, dz=a, nx=node_count, ny=1, nz=1, unit_length=unit_length)
    sim = Sim(mesh, "dw_cobalt")
    sim.Ms = Ms
    sim.set_m(initial_m)
    sim.do_precession = do_precession
    sim.add(UniformExchange(A))
    sim.add(UniaxialAnisotropy(K1, (0, 0, 1)))
    sim.pins = lambda r: 1 if (r[0] < a or r[0] > LENGTH - a) else 0
    return sim


def compute_domain_wall_cobalt(end_time=1e-9):
    sim = setup_domain_wall_cobalt()
    sim.run_until(end_time)
    return sim.mesh.coordinates[:, 0], sim.spin


def test_domain_wall_cobalt():
    xs, m = compute_domain_wall_cobalt()
    m.shape = (-1, 3)
    mz = m[:, 2]
    ref = [reference_mz(x) for x in xs]
    diff = np.max(np.abs(mz - ref))
    print "max difference between simulation and reference: ", diff
    assert diff < 1e-2

if __name__ == '__main__':
    xs, m = compute_domain_wall_cobalt()
    m.shape = (-1, 3)
    mz = m[:, 2]
    ref = [reference_mz(x) for x in xs]
    print "max difference between simulation and reference: ", np.max(np.abs(mz - ref))
    plt.plot(xs, mz, label="finmag")
    plt.plot(xs, ref, label="ref")
    plt.xlabel('x [nm]')
    plt.ylabel('m')
    plt.legend()
    plt.title('domain wall in Co')
    plt.savefig('domain_wall_Co.png')
