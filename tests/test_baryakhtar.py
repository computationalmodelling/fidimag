"""The Landau-Lifshitz-Baryakhtar equation: the kernels against NumPy written from the equations, and the drivers.

Vector fields are x0, y0, z0, x1, ... and the cells are numbered with x running fastest, as everywhere in Fidimag."""
import numpy as np
import pytest

import fidimag.extensions.baryakhtar_clib as clib
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, UniaxialAnisotropy, Zeeman

NX, NY, NZ = 5, 4, 3
N = NX * NY * NZ
DX, DY, DZ = 2e-9, 3e-9, 1.5e-9


def random_state(seed=0, vacuum=True):
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(N, 3))
    m /= np.linalg.norm(m, axis=1)[:, None]
    m *= rng.uniform(0.9, 1.1, size=(N, 1))             # (not unit length: the longitudinal relaxation acts on that)
    h = rng.normal(size=(N, 3)) * 1e5
    Ms = np.full(N, 8e5)
    if vacuum:
        Ms[[3, 17, 18, 40]] = 0.0
    alpha = rng.uniform(0.05, 0.4, size=N)
    pins = np.zeros(N, dtype=np.int32)
    pins[[5, 22]] = 1
    return m, h, Ms, alpha, pins


def laplace(h, Ms):
    """sum over the six neighbours with material of (h_neighbour - h) / d^2, zero in a cell without material"""
    h = h.reshape(NZ, NY, NX, 3)
    has = (Ms > 0).reshape(NZ, NY, NX)
    out = np.zeros_like(h)
    for axis, d in ((2, DX), (1, DY), (0, DZ)):
        for step in (-1, 1):
            neighbour = np.roll(h, -step, axis=axis)
            exists = np.roll(has, -step, axis=axis)
            edge = np.zeros_like(has)
            index = [slice(None)] * 3
            index[axis] = -1 if step == 1 else 0
            edge[tuple(index)] = True                        # (no wrapping: the mesh is not periodic)
            use = has & exists & ~edge
            out += np.where(use[..., None], (neighbour - h) / d ** 2, 0.0)
    return out.reshape(-1)


def test_the_laplacian():
    m, h, Ms, *_ = random_state()
    field = np.zeros(3 * N)
    clib.compute_laplace_field(h.reshape(-1).copy(), field, Ms, DX, DY, DZ, NX, NY, NZ)
    assert np.allclose(field, laplace(h.reshape(-1), Ms), rtol=1e-12, atol=1e-6 * np.abs(field).max() * 1e-6)


def test_the_laplacian_of_a_linear_field_is_zero_inside():
    # a field linear in x has no Laplacian away from the ends of the mesh
    x = (np.arange(N) % NX) * DX
    h = np.zeros((N, 3))
    h[:, 1] = 3.0 * x
    field = np.zeros(3 * N)
    clib.compute_laplace_field(h.reshape(-1).copy(), field, np.full(N, 8e5), DX, DY, DZ, NX, NY, NZ)
    inside = ((np.arange(N) % NX) > 0) & ((np.arange(N) % NX) < NX - 1)
    assert np.abs(field.reshape(-1, 3)[inside]).max() < 1e-3 * 3.0 / DX


def test_the_relaxation_field():
    m, h, Ms, *_ = random_state()
    chi_inv = 1.0 / 2e-3
    field = np.zeros(3 * N)
    clib.compute_relaxation_field(m.reshape(-1).copy(), field, Ms, chi_inv, N)
    mm = np.sum(m ** 2, axis=1)
    assert np.allclose(field.reshape(-1, 3), (Ms * chi_inv / 2.0 * (1 - mm))[:, None] * m, rtol=1e-12, atol=0)


def test_the_perpendicular_field():
    m, h, *_ = random_state()
    perp = np.zeros(3 * N)
    clib.compute_perp_field(m.reshape(-1).copy(), h.reshape(-1).copy(), perp, N)
    mm = np.sum(m ** 2, axis=1)
    expected = h - (np.sum(h * m, axis=1) / mm)[:, None] * m
    assert np.allclose(perp.reshape(-1, 3), expected, rtol=1e-12, atol=1e-7)
    assert np.abs(np.sum(perp.reshape(-1, 3) * m, axis=1)).max() < 1e-7 * np.abs(h).max()    # (it is perpendicular to m)


@pytest.mark.parametrize('precession', [True, False])
def test_the_full_equation(precession):
    m, h, Ms, alpha, pins = random_state()
    delta_h = np.random.default_rng(5).normal(size=(N, 3)) * 1e13
    gamma, beta = 2.211e5, 3e-17
    ydot = np.full(3 * N, np.nan)
    clib.compute_llg_rhs_baryakhtar(ydot, m.reshape(-1).copy(), h.reshape(-1).copy(), delta_h.reshape(-1).copy(), alpha, beta, pins, gamma, N, precession)
    expected = gamma * (alpha[:, None] * h - beta * delta_h)
    if precession:
        expected += -gamma * np.cross(m, h)
    expected[pins > 0] = 0
    assert np.allclose(ydot.reshape(-1, 3), expected, rtol=1e-12, atol=1e-9 * np.abs(expected).max())


@pytest.mark.parametrize('precession', [True, False])
@pytest.mark.parametrize('default_c', [0.0, 5e9, -1.0])
def test_the_reduced_equation(precession, default_c):
    m, h, Ms, alpha, pins = random_state(vacuum=False)
    m /= np.linalg.norm(m, axis=1)[:, None]                  # (unit vectors: the reduced equation assumes that)
    mm = np.sum(m ** 2, axis=1)
    hp = h - np.sum(h * m, axis=1)[:, None] * m
    delta_hp = np.random.default_rng(6).normal(size=(N, 3)) * 1e13
    gamma, beta = 2.211e5, 3e-17
    ydot = np.full(3 * N, np.nan)
    clib.compute_llg_rhs_baryakhtar_reduced(ydot, m.reshape(-1).copy(), hp.reshape(-1).copy(), delta_hp.reshape(-1).copy(), alpha, beta, pins, gamma, N,
                                            precession, default_c)
    # dm/dt = -g m x hp - g m x (m x (alpha hp - beta laplace hp)) [+ c (1 - m.m) m]
    damping = alpha[:, None] * hp - beta * delta_hp
    expected = gamma * (mm[:, None] * damping - np.sum(m * damping, axis=1)[:, None] * m)
    if precession:
        expected += -gamma * np.cross(m, hp)
    c = 6 * np.linalg.norm(expected, axis=1) if default_c < 0 else default_c
    expected = expected + (c * (1 - mm))[:, None] * m
    expected[pins > 0] = 0
    assert np.allclose(ydot.reshape(-1, 3), expected, rtol=1e-12, atol=1e-9 * np.abs(expected).max())


def test_the_reduced_equation_without_a_laplacian_is_the_landau_lifshitz_equation():
    """beta = 0 and the field perpendicular to m: dm/dt = -g m x h - g alpha m x (m x h)"""
    m, h, Ms, alpha, pins = random_state(vacuum=False)
    m /= np.linalg.norm(m, axis=1)[:, None]
    pins[:] = 0
    hp = h - np.sum(h * m, axis=1)[:, None] * m
    gamma = 2.211e5
    ydot = np.zeros(3 * N)
    clib.compute_llg_rhs_baryakhtar_reduced(ydot, m.reshape(-1).copy(), hp.reshape(-1).copy(), np.zeros(3 * N), alpha, 0.0, pins, gamma, N, True, 0.0)
    expected = -gamma * np.cross(m, h) - gamma * alpha[:, None] * np.cross(m, np.cross(m, h))
    assert np.allclose(ydot.reshape(-1, 3), expected, rtol=1e-10, atol=1e-9 * np.abs(expected).max())


# the physics, which does not depend on how the kernels are written ----------------------------------------------------

@pytest.mark.parametrize('beta', [0.0, 1e-17, 1e-15])
def test_the_equation_dissipates_energy(beta):
    """The relaxation terms (alpha h - beta laplace h) move m along +h, the direction in which the energy falls (dE/dt ~ -h . dm/dt), so
    with beta > 0 their power is positive: the discrete Laplacian is negative semi-definite."""
    m, h, Ms, alpha, pins = random_state()
    pins[:] = 0
    h = h.reshape(-1)
    delta_h = np.zeros(3 * N)
    clib.compute_laplace_field(h.copy(), delta_h, Ms, DX, DY, DZ, NX, NY, NZ)
    ydot = np.zeros(3 * N)
    clib.compute_llg_rhs_baryakhtar(ydot, m.reshape(-1).copy(), h.copy(), delta_h, alpha, beta, pins, 2.211e5, N, False)    # (no precession)
    weight = np.repeat(Ms > 0, 3)
    assert np.sum((h * ydot)[weight]) >= 0
    assert np.sum((h * delta_h)[weight]) <= 0       # (the Laplacian is negative semi-definite)


def test_precession_does_no_work():
    m, h, Ms, alpha, pins = random_state(vacuum=False)
    pins[:] = 0
    ydot = np.zeros(3 * N)
    clib.compute_llg_rhs_baryakhtar(ydot, m.reshape(-1).copy(), h.reshape(-1).copy(), np.zeros(3 * N), np.zeros(N), 0.0, pins, 2.211e5, N, True)
    assert np.abs(np.sum(ydot.reshape(-1, 3) * h)) < 1e-9 * np.abs(ydot).max() * np.abs(h).max() * N


# the drivers --------------------------------------------------------------------------------------------------------

def run(driver, directory, monkeypatch, beta=1e-18, chi=1e-3, t=1e-13):
    monkeypatch.chdir(directory)
    mesh = CuboidMesh(nx=NX, ny=NY, nz=NZ, dx=DX * 1e9, dy=DY * 1e9, dz=DZ * 1e9, unit_length=1e-9)
    sim = Sim(mesh, name='llbar', driver=driver)
    sim.set_Ms(lambda r: 0.0 if r[0] > 7e-9 and r[1] < 2e-9 else 8e5)
    sim.driver.gamma = 2.211e5
    sim.alpha = lambda r: 0.05 + 0.05 * r[0] * 1e8             # (alpha that varies from cell to cell)
    sim.driver.beta = beta
    if driver == 'llbar_full':
        sim.driver.chi = chi
    sim.set_m(lambda r: (np.sin(4e8 * r[0]), np.cos(3e8 * r[1]), 0.3 + 1e8 * r[2]))
    sim.add(UniformExchange(1.3e-11))
    sim.add(UniaxialAnisotropy(4e4, axis=(0, 0, 1)))
    sim.add(Zeeman((1e4, 5e3, 0)))
    sim.driver.set_tols(rtol=1e-8, atol=1e-10)
    sim.driver.run_until(t)
    return sim.spin.copy()


@pytest.mark.parametrize('driver', ['llbar', 'llbar_full'])
def test_the_driver_runs_and_is_reproducible(driver, tmp_path, monkeypatch):
    """(`llbar` could not be made, and the kernels read past the end of alpha, so that two runs differed)"""
    runs = []
    for i in range(2):
        (tmp_path / str(i)).mkdir()
        runs.append(run(driver, tmp_path / str(i), monkeypatch))
    assert np.all(np.isfinite(runs[0])) and np.array_equal(runs[0], runs[1])


def test_the_laplacian_term_changes_the_dynamics(tmp_path, monkeypatch):
    (tmp_path / 'a').mkdir(), (tmp_path / 'b').mkdir()
    without = run('llbar', tmp_path / 'a', monkeypatch, beta=0.0)
    with_beta = run('llbar', tmp_path / 'b', monkeypatch, beta=1e-15)
    assert np.abs(without - with_beta).max() > 1e-8


def test_a_uniform_length_relaxes_at_the_rate_alpha_gamma_Ms_over_chi(tmp_path, monkeypatch):
    """With no field but the relaxation one, dm/dt = gamma alpha (Ms / 2 chi) (1 - m.m) m, so a length 1 - eps decays as exp(-gamma alpha Ms t / chi)."""
    monkeypatch.chdir(tmp_path)
    mesh = CuboidMesh(nx=2, ny=2, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh, name='length', driver='llbar_full')
    Ms, alpha, gamma, chi, eps = 8e5, 0.1, 2.211e5, 1e-3, 1e-2
    sim.set_Ms(Ms)
    sim.driver.gamma, sim.alpha, sim.driver.beta, sim.driver.chi = gamma, alpha, 0.0, chi
    sim.set_m((0, 0, 1 - eps), normalise=False)
    sim.driver.set_tols(rtol=1e-10, atol=1e-14)
    rate = gamma * alpha * Ms / chi
    t = 3.0 / rate
    sim.driver.run_until(t)
    length = np.linalg.norm(sim.spin.reshape(-1, 3), axis=1)
    assert np.allclose(1 - length, eps * np.exp(-rate * t), rtol=2e-2)
