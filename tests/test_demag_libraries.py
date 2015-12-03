"""
Testing of the different functions to compute the demag
field for atomistic simulations (in micromagnetic simulations
we only use cuboid meshes, so we only have the FFT approach
for the dipolar field)
"""

from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import Demag, DemagFull, DemagHexagonal
from fidimag.atomistic import Constant

const = Constant()


def test_cuboid_demags():
    """
    Test a brute force calculation of the demagnetising field, called
    DemagFull, based on the sum of the dipolar contributions of the whole
    system for every lattice site, against the default FFT approach for the
    demag field. We compute the energies scaled in meV.
    This test is performed in a cuboid mesh to assure that the DemagFull
    library is calculating the same than the default demag function
    """

    mesh = CuboidMesh(0.27, 0.27, 0.27, 4, 4, 1, unit_length=1e-9)
    mu_s = 2 * const.mu_B

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m((0, 0, 1))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('demag_full').compute_field()
    # print sim.get_interaction('demag_full').field
    demag_full_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m((0, 0, 1))

    sim2.add(Demag())
    sim2.get_interaction('demag').compute_field()
    sim2.compute_energy()

    demag_fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (demag_full_energy - demag_fft_energy) < 1e-10


def test_hexagonal_demags():
    """
    Comparison of the FFT approach for hexagonal meshes, named
    DemagHexagonal, where it is used a system with the double number
    of nodes along the x direction (i.e. a mesh with twice the number
    of nodes of the original mesh), against the full calculation
    of the Demag field
    """

    # It seems that using the square and diagonal
    # alignment of the spins positions in the lattice
    # produce the same result, although the field is
    # diferent. The library is based on the square alignment
    mesh = HexagonalMesh(0.27 * 0.5, 4, 4, 
                         unit_length=1e-9, 
                         alignment='square')
    mu_s = 2 * const.mu_B

    sim = Sim(mesh)
    sim.mu_s = mu_s

    sim.set_m((0, 0, 1))
    # Brute force demag calculation
    sim.add(DemagFull())

    sim.get_interaction('demag_full').compute_field()
    print sim.get_interaction('demag_full').field
    demag_full_energy = sim.compute_energy() / const.meV

    # Demag using the FFT approach and a larger mesh
    sim2 = Sim(mesh)
    sim2.mu_s = mu_s

    sim2.set_m((0, 0, 1))

    sim2.add(DemagHexagonal())
    sim2.get_interaction('demag_hex').compute_field()
    sim2.compute_energy()

    demag_2fft_energy = sim2.compute_energy() / const.meV

    # We compare both energies scaled in meV
    assert (demag_full_energy - demag_2fft_energy) < 1e-10


if __name__ == '__main__':
    test_cuboid_demags()
    test_hexagonal_demags()
