from fidimag.atomistic import Anisotropy, CubicAnisotropy
from fidimag.common import CuboidMesh
from fidimag.atomistic import Sim
import numpy as np


def test_anis():
    mesh = CuboidMesh(nx=5, ny=3, nz=2)
    spin = np.zeros(90)
    anis = Anisotropy(Ku=1, axis=[1, 0, 0])
    # Magnetisation is one value per mesh node
    mu_s = np.ones(30)
    mu_s_inv = np.ones(30)
    anis.setup(mesh, spin, mu_s, mu_s_inv)
    field = anis.compute_field()
    assert len(mesh.coordinates) == 5 * 3 * 2
    assert np.max(field) == 0
    spin[0] = 99
    field = anis.compute_field()
    assert field[0] == 2 * 99

def test_anis_cubic_macrospin():
    """Test 1D system anisotropy
    """
    mesh = CuboidMesh(nx=1, ny=1, nz=1)
    spin = np.array([0.6,0.8,0])
    anis = CubicAnisotropy(Kc=1.23)
    # Magnetisation is one value per mesh node
    mu_s = np.ones(1)
    mu_s_inv = np.ones(1)
    anis.setup(mesh, spin, mu_s, mu_s_inv)
    field = anis.compute_field()
    assert  np.max(field-np.array([1.06272, 2.51904, 0.]))<1e-6
    # OLD CODE:
    # assert  np.max(field-np.array([-1.06272, -2.51904, 0.]))<1e-6
    # assert abs(energy-0.663216)<1e-5
    energy = anis.compute_energy()
    assert abs(energy - (-0.663216)) < 1e-5

def test_anis_cubic():
    """Test 2D system anisotropy

                                    x, y, z    x, y, z
    Check only two spins at sites: (0, 0, 0), (0, 1, 0)
    """
    mesh = CuboidMesh(nx=2, ny=2, nz=1)
    # Do not normalise for the test:
    spin = np.array([0.6, 0.8, 0,
                     0.3, 0.3, 0.3,
                     0.5, 0.5, 0.2,
                     0.3, 0.3, 0.3,
                     ])
    anis = CubicAnisotropy(Kc=np.array([1.23, 0.88, 0.88, 1.23]))
    # Magnetisation is one value per mesh node
    mu_s = np.ones(4)
    mu_s_inv = np.ones(4)
    anis.setup(mesh, spin, mu_s, mu_s_inv)
    field = (anis.compute_field()).reshape(-1, 3)
    assert  np.max(field[0] - np.array([1.06272, 2.51904, 0.])) < 1e-6
    assert  np.max(field[2] - np.array([0.44, 0.44, 0.02816])) < 1e-6
    energy = anis.compute_energy()
    assert abs(anis.energy[0] - (-0.663216)) < 1e-5
    assert abs(anis.energy[2] - (-0.111408)) < 1e-5
    # print(anis.energy)


if __name__ == '__main__':
    test_anis()
    test_anis_cubic()
    test_anis_cubic_macrospin()
