"""
Testing of the different functions to compute the demag
field for atomistic simulations (in micromagnetic simulations
we only use cuboid meshes, so we only have the FFT approach
for the dipolar field)
"""
import fidimag
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import Demag, DemagFull, DemagHexagonal
import fidimag.common.constant as const
import numpy as np


def m_init_dw(pos, N, a):
    """
    This creates a pseudo domain wall along x
    N is the number of atoms and a the lattice spacing
    """
    x = 2 * pos[0] / (N * a) - 1
    mx = x
    mz = np.sqrt(1 - x ** 2.)
    return (mx, 0, mz)


def m_init_2Dvortex(pos, centre):
    x, y = pos[0] - centre[0], pos[1] - centre[1]

    if np.sqrt(x ** 2 + y ** 2) <= 2.5:
        # Polar coordinates:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        # This determines the profile we want for the
        # skyrmion-like vortex
        # Single twisting: k = pi / R
        k = np.pi / 2.5

        # We define here a 'hedgehog' skyrmion pointing up
        return (np.sin(k * r) * np.cos(phi),
                np.sin(k * r) * np.sin(phi),
                np.cos(k * r))
    else:
        return (0, 0, -1)


def test_hexagonal_demags_1Dchain():
    """
    Comparison of the FFT approach for hexagonal meshes, named
    DemagHexagonal, where it is used a system with the double number
    of nodes along the x direction (i.e. a mesh with twice the number
    of nodes of the original mesh), against the full calculation
    of the Demag field
    """
    # Number of atoms
    N = 12
    a = 0.4
    mesh = HexagonalMesh(a * 0.5, N, 1,
                         unit_length=1e-9,
                         alignment='square')
    mu_s = 2 * const.mu_B

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m(lambda pos: m_init_dw(pos, N, a))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('DemagFull').compute_field()
    sim.get_interaction('DemagFull').field
    DemagFull_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach and a larger mesh
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m(lambda pos: m_init_dw(pos, N, a))

    sim2.add(DemagHexagonal())
    sim2.get_interaction('DemagHexagonal').compute_field()
    sim2.compute_energy()

    demag_2fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (DemagFull_energy - demag_2fft_energy) < 1e-10


def test_cuboid_demags_1Dchain():
    """
    Test a brute force calculation of the demagnetising field, called
    DemagFull, based on the sum of the dipolar contributions of the whole
    system for every lattice site, against the default FFT approach for the
    demag field. We compute the energies scaled in meV.
    This test is performed in a cuboid mesh to assure that the DemagFull
    library is calculating the same than the default demag function
    """
    N = 12
    a = 0.4
    mesh = CuboidMesh(a, a, a, N, 1, 1, unit_length=1e-9)
    mu_s = 2 * const.mu_B

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m(lambda pos: m_init_dw(pos, N, a))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('DemagFull').compute_field()
    # print sim.get_interaction('DemagFull').field
    DemagFull_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m(lambda pos: m_init_dw(pos, N, a))

    sim2.add(Demag())
    sim2.get_interaction('Demag').compute_field()
    sim2.compute_energy()

    demag_fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (DemagFull_energy - demag_fft_energy) < 1e-10


def test_cuboid_demags_2D():
    """
    Comparison of the FFT approach for hexagonal meshes, named
    DemagHexagonal, where it is used a system with the double number
    of nodes along the x direction (i.e. a mesh with twice the number
    of nodes of the original mesh), against the full calculation
    of the Demag field
    """
    # Number of atoms
    N = 15
    a = 0.4

    mesh = CuboidMesh(a, a, a, N, N, 1, unit_length=1e-9)
    mu_s = 2 * const.mu_B

    # Centre
    xc = (mesh.Lx * 0.5)
    yc = (mesh.Ly * 0.5)

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m(lambda pos: m_init_2Dvortex(pos, (xc, yc)))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('DemagFull').compute_field()
    # print sim.get_interaction('DemagFull').field
    DemagFull_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m(lambda pos: m_init_2Dvortex(pos, (xc, yc)))

    sim2.add(Demag())
    sim2.get_interaction('Demag').compute_field()
    sim2.compute_energy()

    demag_fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (DemagFull_energy - demag_fft_energy) < 1e-10


def test_hexagonal_demags_2D():
    """
    Comparison of the FFT approach for hexagonal meshes, named
    DemagHexagonal, where it is used a system with the double number
    of nodes along the x direction (i.e. a mesh with twice the number
    of nodes of the original mesh), against the full calculation
    of the Demag field
    """
    # Number of atoms
    N = 15
    a = 0.4
    mesh = HexagonalMesh(a * 0.5, N, N,
                         unit_length=1e-9,
                         alignment='square')

    # Centre
    xc = (mesh.Lx * 0.5)
    yc = (mesh.Ly * 0.5)

    mu_s = 2 * const.mu_B

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m(lambda pos: m_init_2Dvortex(pos, (xc, yc)))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('DemagFull').compute_field()
    sim.get_interaction('DemagFull').field
    DemagFull_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach and a larger mesh
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m(lambda pos: m_init_2Dvortex(pos, (xc, yc)))

    sim2.add(DemagHexagonal())
    sim2.get_interaction('DemagHexagonal').compute_field()
    sim2.compute_energy()

    demag_2fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (DemagFull_energy - demag_2fft_energy) < 1e-10



def test_demag_2d_pbc():
    """
    Attempt to check that demag with 2d_pbc option does
    not give a nonsensical answer.
    """
    A=1.3e-11
    Ms=8.6e5
    n = 40
    d = 2.5
    
    mesh = fidimag.common.CuboidMesh(nx=n, ny=n, nz=1, dx=d, dy=d, dz=d, unit_length=1e-9, periodicity=(True, True, False))
    sim = fidimag.micro.Sim(mesh, name="pbc_2d_bug")
    sim.set_Ms(Ms)
    sim.set_m((0, 0, 1.0), normalise=True)
    demag = fidimag.micro.Demag(pbc_2d=True)

    sim.add(demag)
    sim.compute_effective_field(0)

    assert not np.isnan(demag.demag.tensor_xx).any()
    assert not np.isnan(demag.demag.tensor_xy).any()
    assert not np.isnan(demag.demag.tensor_xz).any()
    assert not np.isnan(demag.demag.tensor_yy).any()
    assert not np.isnan(demag.demag.tensor_yz).any()
    assert not np.isnan(demag.demag.tensor_zz).any()
    assert not np.isnan(sim.field).any(), "NaN in demag array"

if __name__ == '__main__':
    test_hexagonal_demags_1Dchain()
    test_cuboid_demags_1Dchain()
    test_hexagonal_demags_2D()
    test_cuboid_demags_2D()
