"""The classical mean-field driver, checked against results that can be derived by hand."""
import numpy as np
import pytest

import fidimag.common.constant as const
import fidimag.extensions.clib as clib
from fidimag.atomistic import Sim, Exchange, DMI, Zeeman, UniformExchange, Anisotropy
from fidimag.common import CuboidMesh
from fidimag.atomistic.cmf import CMFResult

kB = const.k_B
J = 100 * kB                 # for a simple cubic lattice (z = 6) T_c = z J / (3 kB) = 200 K
TC = 6 * J / (3 * kB)
MU = const.mu_B


def langevin(x):
    return 1 / np.tanh(x) - 1 / x


def ferromagnet(n=4, mu_s=MU):
    mesh = CuboidMesh(nx=n, ny=n, nz=n, periodicity=(True, True, True))
    sim = Sim(mesh, name='cmf_ferro', driver='cmf')
    sim.set_mu_s(mu_s)
    sim.add(Exchange(J))
    return sim


def magnetisation(sim):
    return np.linalg.norm(sim.spin.reshape(-1, 3).mean(axis=0))


def target(field_z, T, mu_s=1.0, pinned=0):
    spin, resid = np.zeros(3), np.zeros(3)
    field = np.array([0.0, 0.0, field_z])
    clib.compute_cmf_residual(spin, field, np.array([mu_s]), np.array([pinned], dtype=np.int32), resid, T, 1.0)
    return resid.copy()


# ---------------------------------------------------------------------------- the kernel


def test_the_langevin_function_has_its_known_values():
    assert target(1.0, 1.0)[2] == pytest.approx(1 / np.tanh(1) - 1, rel=1e-13)       # L(1) = 0.31303528...
    assert target(1.0, 1.0)[2] == pytest.approx(0.3130352854993, rel=1e-12)
    assert target(100.0, 1.0)[2] == pytest.approx(1 - 1 / 100, rel=1e-13)            # L(x) -> 1 - 1/x


def test_the_response_to_a_vanishing_field_is_one_third_and_is_continuous_at_the_series_switch():
    assert target(1e-6, 1.0)[2] / 1e-6 == pytest.approx(1 / 3, rel=1e-9)
    a, b = target(0.1 - 1e-12, 1.0)[2], target(0.1 + 1e-12, 1.0)[2]
    assert a == pytest.approx(b, rel=1e-9)
    # the series is L(x)/x = 1/3 - x^2/45 + 2 x^4/945
    x = 0.05
    assert target(x, 1.0)[2] / x == pytest.approx(1 / 3 - x ** 2 / 45 + 2 * x ** 4 / 945, rel=1e-9)


def test_the_moment_is_along_the_field_and_shorter_than_one():
    spin, resid = np.zeros(3), np.zeros(3)
    field = np.array([0.3, -0.2, 0.5])
    clib.compute_cmf_residual(spin, field, np.array([1.0]), np.zeros(1, dtype=np.int32), resid, 1.0, 1.0)
    assert np.allclose(np.cross(resid, field), 0, atol=1e-14) and resid @ field > 0
    assert 0 < np.linalg.norm(resid) < 1
    assert np.linalg.norm(resid) == pytest.approx(langevin(np.linalg.norm(field)), rel=1e-12)


def test_at_zero_temperature_the_moment_is_a_unit_vector_along_the_field():
    assert np.allclose(target(0.3, 0.0), [0, 0, 1])


def test_pinned_spins_and_spins_without_a_moment_are_left_alone():
    assert np.all(target(1.0, 1.0, pinned=1) == 0)
    assert np.all(target(1.0, 1.0, mu_s=0.0) == 0)


# ---------------------------------------------------------------------------- the ferromagnet


def test_above_the_critical_temperature_the_moment_vanishes():
    sim = ferromagnet()
    sim.set_m((0.3, 0.1, 1.0))
    result = sim.driver.minimise(T=1.2 * TC, history=8, stopping_dm=1e-10)
    assert result.converged and magnetisation(sim) < 1e-8


def test_the_susceptibility_follows_the_curie_weiss_law():
    # M = mu B / (3 kB (T - T_c)) for a field small enough to be linear: the field is mu B / kB T = 2e-3 here
    T = 2 * TC
    B = 2e-3 * kB * T / MU
    sim = ferromagnet()
    sim.add(Zeeman((0, 0, B)))
    sim.set_m((0, 0, 1))
    assert sim.driver.minimise(T=T, history=8, stopping_dm=1e-12).converged
    assert magnetisation(sim) == pytest.approx(MU * B / (3 * kB * (T - TC)), rel=1e-4)


def test_below_the_critical_temperature_the_moment_solves_the_mean_field_equation():
    # M = L(3 M T_c / T)
    sim = ferromagnet()
    for t in (0.3, 0.6, 0.9, 0.99):
        sim.set_m((0.3, 0.1, 1.0))
        assert sim.driver.minimise(T=t * TC, history=8, stopping_dm=1e-12).converged
        M = magnetisation(sim)
        assert M == pytest.approx(langevin(3 * M / t), rel=1e-8)
        # every site holds the same moment
        assert np.ptp(np.linalg.norm(sim.spin.reshape(-1, 3), axis=1)) < 1e-9


def test_close_to_the_critical_temperature_the_order_parameter_has_the_square_root_law():
    # L(x) = x/3 - x^3/45 gives M^2 = (5/3) (1 - T/T_c) to leading order
    sim = ferromagnet()
    sim.set_m((0, 0, 1))
    t = 0.999
    assert sim.driver.minimise(T=t * TC, history=10, stopping_dm=1e-12).converged
    assert magnetisation(sim) ** 2 == pytest.approx(5 / 3 * (1 - t), rel=0.01)


def test_at_zero_temperature_the_spins_align_and_the_energy_is_the_exchange_one():
    sim = ferromagnet()
    sim.set_m(lambda r: (0.1 * np.sin(r[0] * 3e9), 0.1, 1.0))
    result = sim.driver.minimise(T=0, history=8, stopping_dm=1e-12)
    assert result.converged
    assert np.allclose(np.linalg.norm(sim.spin.reshape(-1, 3), axis=1), 1)
    assert result.total_energy / sim.driver.n == pytest.approx(-3 * J, rel=1e-9)     # -z J / 2 per site


def test_the_moment_does_not_depend_on_the_direction_it_starts_in():
    sim = ferromagnet()
    values = []
    for direction in ((0, 0, 1), (1, 0, 0), (1, 1, 1)):
        sim.set_m(direction)
        assert sim.driver.minimise(T=0.8 * TC, history=8, stopping_dm=1e-12).converged
        values.append(magnetisation(sim))
    assert np.ptp(values) < 1e-8


def test_anderson_mixing_finds_the_same_state_in_fewer_steps():
    sim = ferromagnet()
    sim.set_m((0, 0, 1))
    plain = sim.driver.minimise(T=0.99 * TC, stopping_dm=1e-10, mixing=0.5)
    m_plain = magnetisation(sim)
    sim.set_m((0, 0, 1))
    anderson = sim.driver.minimise(T=0.99 * TC, stopping_dm=1e-10, mixing=0.5, history=8)
    assert plain.converged and anderson.converged
    assert magnetisation(sim) == pytest.approx(m_plain, abs=1e-7)
    assert anderson.n_iterations * 5 < plain.n_iterations


# ---------------------------------------------------------------------------- DMI


def dmi_chain(D_over_J, n=64):
    mesh = CuboidMesh(nx=n, ny=1, nz=1, periodicity=(True, False, False))
    sim = Sim(mesh, name='cmf_dmi', driver='cmf')
    sim.set_mu_s(MU)
    sim.add(Exchange(J))
    sim.add(DMI(D_over_J * J, dmi_type='bulk'))
    return sim


def test_a_chain_with_dmi_orders_below_and_not_above_the_temperature_where_the_paramagnet_becomes_unstable():
    # The paramagnet is unstable when 3 kB T = lambda_max, with lambda_max = 2 sqrt(J^2 + D^2) for a chain
    # (the largest eigenvalue of the coupling, at tan q = D / J)
    for d in (0.3, 1.0):
        T_instability = 2 * np.sqrt(J ** 2 + (d * J) ** 2) / (3 * kB)
        sim = dmi_chain(d)
        for factor, ordered in ((0.97, True), (1.03, False)):
            sim.set_m(np.random.default_rng(1).normal(size=3 * 64) * 0.05)
            sim.driver.minimise(T=factor * T_instability, history=10, stopping_dm=1e-6, max_steps=3000)
            assert (np.linalg.norm(sim.spin.reshape(-1, 3), axis=1).max() > 0.05) == ordered


def test_the_helix_has_the_pitch_tan_q_equals_d_over_j_and_the_energy_minus_the_root_of_j_squared_plus_d_squared():
    # D = J gives q = 2 pi 8 / 64 on 64 sites; the energy per bond is -sqrt(J^2 + D^2) = -sqrt(2) J
    n = 64
    sim = dmi_chain(1.0, n)
    x = np.arange(n)
    for k in (8, -8):
        q = 2 * np.pi * k / n
        m = np.stack([np.zeros(n), np.sin(q * x), np.cos(q * x)], axis=1) + 0.2 * np.random.default_rng(0).normal(size=(n, 3))
        sim.set_m((m / np.linalg.norm(m, axis=1)[:, None]).ravel())
        result = sim.driver.minimise(T=1e-3 * J / kB, history=10, stopping_dm=1e-10)
        assert result.converged
        assert np.argmax(np.abs(np.fft.fft(sim.spin.reshape(-1, 3)[:, 1]))[1:n // 2]) + 1 == 8
        # (the thermal shortening of the moment, |S|^2, takes 1e-3 off)
        mean_s2 = np.mean(np.sum(sim.spin.reshape(-1, 3) ** 2, axis=1))
        assert result.total_energy / n / J == pytest.approx(-np.sqrt(2) * mean_s2, rel=1e-6)


# ---------------------------------------------------------------------------- the driver


def test_single_ion_anisotropy_is_refused():
    sim = ferromagnet()
    sim.add(Anisotropy(1e-24, axis=(0, 0, 1)))
    with pytest.raises(NotImplementedError, match="anisotropy"):
        sim.driver.minimise(T=100)


def test_a_run_that_does_not_converge_says_so_and_does_not_raise(caplog):
    sim = ferromagnet()
    sim.set_m((0, 0, 1))
    result = sim.driver.minimise(T=0.999 * TC, max_steps=5, stopping_dm=1e-14)
    assert isinstance(result, CMFResult)
    assert not result.converged and result.n_iterations == 5 and "max_steps" in result.reason


def test_pinned_spins_do_not_move():
    sim = ferromagnet()
    sim.set_m(lambda r: (0.3, 0, 1) if r[0] < 1.5 else (1, 0, 0))
    sim.set_pins(lambda r: 1 if r[0] < 1.5 else 0)
    before = sim.spin.copy()
    pinned = np.repeat(sim._pins > 0, 3)
    assert 0 < pinned.sum() < pinned.size
    sim.driver.minimise(T=0.5 * TC, history=4, max_steps=50)
    assert np.array_equal(sim.spin[pinned], before[pinned])


@pytest.mark.parametrize('kwargs', [dict(T=-1), dict(T=1, mixing=0), dict(T=1, mixing=1.5), dict(T=1, history=33)])
def test_wrong_arguments_are_refused(kwargs):
    with pytest.raises(ValueError):
        ferromagnet().driver.minimise(**kwargs)
