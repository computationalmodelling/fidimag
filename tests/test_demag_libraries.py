"""
Testing a brute force calculation of the demagnetising field
against the default approach, based on a Fourier transformation.
"""

from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import Demag, DemagFull
from fidimag.atomistic import Constant

const = Constant()


def test_cuboid_demags():

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

if __name__ == '__main__':
    test_cuboid_demags()
