"""
Every interaction reports the same thing in `energy` and `total_energy`.

`energy[i]` is the energy of cell (or site) `i` in joules, and
`total_energy` is their sum. Before this was pinned down the array meant an
energy density in most micromagnetic classes, an already weighted energy in
`micro/demag.py`, and nothing at all in `micro/zeeman.py`, which never wrote
to it; two of the atomistic demagnetising classes returned a total without
ever storing it. Anything summing or differencing per cell energies across
interactions was silently wrong.
"""
import numpy as np
import pytest

import fidimag
from fidimag.common import CuboidMesh
from fidimag.common.constant import k_B, mu_B


Ms = 8.6e5
A = 1.3e-11
Ku = 1e5
D = 1e-3


def _micro_sim():
    mesh = CuboidMesh(nx=6, ny=5, nz=2, dx=3, dy=3, dz=3, unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, name='per_cell_micro')
    sim.set_Ms(Ms)
    np.random.seed(3)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    return sim


def _atom_sim():
    mesh = CuboidMesh(nx=6, ny=5, nz=2)
    sim = fidimag.atomistic.Sim(mesh, name='per_cell_atom')
    sim.set_mu_s(2 * mu_B)
    np.random.seed(3)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    return sim


MICRO = [
    lambda: fidimag.micro.UniformExchange(A),
    lambda: fidimag.micro.Demag(),
    lambda: fidimag.micro.SimpleDemag(),
    lambda: fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)),
    lambda: fidimag.micro.UniaxialAnisotropy4(K1=Ku, K2=0.1 * Ku,
                                              axis=(0, 0, 1)),
    lambda: fidimag.micro.ExchangeRKKY(sigma=1e-4),
    lambda: fidimag.micro.DMI(D=D, dmi_type='bulk'),
    lambda: fidimag.micro.DMI(D=D, dmi_type='interfacial'),
    lambda: fidimag.micro.Zeeman((0, 0, 1e4)),
]

ATOM = [
    lambda: fidimag.atomistic.UniformExchange(J=50 * k_B),
    lambda: fidimag.atomistic.Anisotropy(5 * k_B, axis=(0, 0, 1)),
    lambda: fidimag.atomistic.CubicAnisotropy(Kc=5 * k_B),
    lambda: fidimag.atomistic.DMI(D=1 * k_B),
    lambda: fidimag.atomistic.Zeeman((0, 0, 0.1)),
    lambda: fidimag.atomistic.Demag(),
    lambda: fidimag.atomistic.DemagFull(),
    lambda: fidimag.atomistic.DemagFMM(order=4, ncrit=64, theta=0.5),
]


def _check(sim, make):
    interaction = make()
    sim.add(interaction)
    total = interaction.compute_energy()

    assert isinstance(total, float) or np.isscalar(total), (
        '%s returned %r rather than a total'
        % (interaction.name, type(total)))
    assert interaction.energy.shape == (sim.mesh.n,)
    assert np.all(np.isfinite(interaction.energy))
    assert interaction.total_energy == pytest.approx(total, rel=1e-12)
    assert np.sum(interaction.energy) == pytest.approx(total, rel=1e-10,
                                                       abs=1e-30)
    return total


@pytest.mark.parametrize('make', MICRO,
                         ids=lambda f: f().__class__.__name__)
def test_micro_energy_sums_to_the_total(make):
    _check(_micro_sim(), make)


@pytest.mark.parametrize('make', ATOM,
                         ids=lambda f: f().__class__.__name__)
def test_atomistic_energy_sums_to_the_total(make):
    _check(_atom_sim(), make)


def test_driver_total_matches_the_sum_over_cells():
    """The whole point: per cell energies can be added across interactions."""
    sim = _micro_sim()
    sim.add(fidimag.micro.UniformExchange(A))
    sim.add(fidimag.micro.Demag())
    sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))
    sim.add(fidimag.micro.Zeeman((0, 0, 1e4)))

    per_cell = np.zeros(sim.mesh.n)
    total = 0.0
    for obj in sim.driver.interactions:
        total += obj.compute_energy()
        per_cell += obj.energy

    assert np.sum(per_cell) == pytest.approx(total, rel=1e-10)
    assert sim.compute_energy() == pytest.approx(total, rel=1e-10)


def test_hexagonal_demag_energy_sums_to_the_total():
    """The hexagonal demagnetising class also returned a total it never kept."""
    from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
    mesh = HexagonalMesh(1, 4, 4)
    sim = fidimag.atomistic.Sim(mesh, name='per_cell_hex')
    sim.set_mu_s(2 * mu_B)
    np.random.seed(3)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    demag = fidimag.atomistic.DemagHexagonal()
    sim.add(demag)

    total = demag.compute_energy()
    assert np.isscalar(total)
    assert demag.total_energy == pytest.approx(total, rel=1e-12)
    assert np.sum(demag.energy) == pytest.approx(total, rel=1e-10, abs=1e-30)


def test_field_is_current_after_compute_energy():
    """
    `compute_energy` leaves `field` at the current spins, not the previous.

    The drivers assemble the effective field by calling `compute_energy` on
    each interaction and then reading its `field`, so anything that reports an
    energy without refreshing its field contributes the field of whatever
    configuration came before. Both atomistic demagnetising classes did that,
    which left the atomistic minimiser relaxing against a demagnetising field
    it never recomputed.
    """
    mesh = CuboidMesh(nx=8, ny=6, nz=1, dx=0.5, dy=0.5, dz=0.5,
                      unit_length=1e-9)
    sim = fidimag.atomistic.Sim(mesh, name='field_current')
    sim.set_mu_s(2 * mu_B)
    np.random.seed(1)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    demag = fidimag.atomistic.Demag()
    sim.add(demag)

    demag.compute_field()
    stale = demag.field.copy()
    assert np.abs(stale).max() > 1e-3, 'demag too weak to detect a lag'

    np.random.seed(9)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    demag.compute_energy()
    after = demag.field.copy()
    truth = demag.compute_field().copy()

    assert not np.allclose(stale, truth), 'the two states must differ'
    assert np.allclose(after, truth)


def test_driver_effective_field_is_right_on_the_first_call():
    """The same thing, through the path the minimisers actually use."""
    mesh = CuboidMesh(nx=8, ny=6, nz=1, dx=0.5, dy=0.5, dz=0.5,
                      unit_length=1e-9)
    sim = fidimag.atomistic.Sim(mesh, name='drv_field',
                                driver='hubert_minimiser')
    sim.set_mu_s(2 * mu_B)
    np.random.seed(1)
    sim.set_m(lambda pos: tuple(np.random.rand(3) - 0.5))
    sim.add(fidimag.atomistic.UniformExchange(J=50 * k_B))
    sim.add(fidimag.atomistic.Demag())

    driver = sim.driver
    driver.compute_effective_field()
    assembled = driver.field.copy()

    fresh = np.zeros_like(assembled)
    for obj in driver.interactions:
        fresh += obj.compute_field()

    assert np.allclose(assembled, fresh)
