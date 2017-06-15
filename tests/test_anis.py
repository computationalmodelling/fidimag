from fidimag.atomistic import Anisotropy, CubicAnisotropy
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
import numpy as np


def test_anis():
    mesh = CuboidMesh(nx=5, ny=3, nz=2)
    spin = np.zeros(90)
    anis = Anisotropy(Ku=1, axis=[1, 0, 0])
    mu_s = np.ones(90)
    anis.setup(mesh, spin, mu_s)
    field = anis.compute_field()
    assert len(mesh.coordinates) == 5 * 3 * 2
    assert np.max(field) == 0
    spin[0] = 99
    field = anis.compute_field()
    assert field[0] == 2 * 99

def test_anis_cubic():
    mesh = CuboidMesh(nx=1, ny=1, nz=1)
    spin = np.array([0.6,0.8,0])
    anis = CubicAnisotropy(Kc=1.23)
    mu_s = np.ones(3)
    anis.setup(mesh, spin, mu_s)
    field = anis.compute_field()
    print field
    assert  np.max(field-np.array([-1.06272,-2.51904, -0.]))<1e-6
    energy = anis.compute_energy()
    assert abs(energy-0.663216)<1e-5


if __name__ == '__main__':
    test_anis()
    test_anis_cubic()
