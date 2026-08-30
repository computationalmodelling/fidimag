"""The analytical Jacobian must be the derivative of the right hand side.

Issue #21: `use_jac=True` appeared to hang. It did not: the Jacobian-times-
vector product was wrong, so GMRES never converged and spent its whole
iteration budget on every solve, half a million evaluations for a single
picosecond of simulated time. Nothing detected that, because nothing compared
the analytical product against the derivative it claims to be.

These tests do compare it, against a finite difference of the right hand side,
which is the only reference that cannot itself be wrong in the same way.
"""
import warnings

import numpy as np
import pytest

from fidimag.common import CuboidMesh
from fidimag.common.constant import mu_B, k_B
from fidimag.micro import Sim as MicroSim, UniformExchange, Demag, Zeeman
from fidimag.atomistic import Sim as AtomSim, \
    UniformExchange as AtomExchange, Anisotropy as AtomAnisotropy


def _random_unit_field(n, seed):
    rng = np.random.RandomState(seed)
    m = rng.normal(size=(n, 3))
    m /= np.linalg.norm(m, axis=1)[:, None]
    return m.reshape(-1)


def _relative_error(sim, seed=7):
    """Compare the analytical J*v with (f(m + eps v) - f(m)) / eps."""
    driver = sim.driver
    m = sim.spin.copy()
    n = m.size

    rng = np.random.RandomState(seed)
    v = rng.normal(size=n)
    v /= np.linalg.norm(v)

    def rhs_at(state):
        # sundials_rhs reads self.spin rather than its argument: in the real
        # integrator the two share memory
        out = np.zeros(n)
        sim.spin[:] = state
        driver.sundials_rhs(0.0, sim.spin, out)
        return out

    fy = rhs_at(m)

    Jv = np.zeros(n)
    sim.spin[:] = m
    driver.sundials_jtimes(v.copy(), Jv, 0.0, m.copy(), fy.copy())

    # the step that balances truncation against round off for this problem
    finite_difference = (rhs_at(m + 1e-7 * v) - fy) / 1e-7
    sim.spin[:] = m

    return (np.linalg.norm(Jv - finite_difference)
            / np.linalg.norm(finite_difference))


def _micro_sim(with_demag=False, with_zeeman=False, do_precession=True):
    mesh = CuboidMesh(3, 3, 10 / 3.0, 4, 4, 4, unit_length=1e-9)
    sim = MicroSim(mesh, name='jac_micro', integrator='cvode_bdf',
                   driver='llg', use_jac=True)
    sim.Ms = 0.86e6
    sim.driver.alpha = 0.5
    sim.driver.do_precession = do_precession
    sim.set_m(_random_unit_field(mesh.n, seed=42))
    sim.add(UniformExchange(A=13e-12))
    if with_demag:
        sim.add(Demag())
    if with_zeeman:
        sim.add(Zeeman([0, 0, 0.1 / (4 * np.pi * 1e-7)]))
    return sim


@pytest.mark.parametrize('do_precession', [False, True])
def test_micromagnetic_jacobian_matches_finite_differences(do_precession):
    sim = _micro_sim(do_precession=do_precession)
    error = _relative_error(sim)
    print('micro precession={}: {:.2e}'.format(do_precession, error))
    assert error < 1e-5


def test_jacobian_with_a_constant_field():
    """A Zeeman field does not vary with m and must not enter the Jacobian.

    That is what the `jac` flag on each interaction is for, and getting it
    wrong would show up here as a term that finite differences do not see.
    """
    sim = _micro_sim(with_zeeman=True)
    assert _relative_error(sim) < 1e-5


def test_jacobian_with_a_non_uniform_damping():
    """alpha is one value per site, and was being indexed per component.

    With a uniform alpha the mistake is invisible; with a varying one it is
    not, and it also read past the end of the array.
    """
    sim = _micro_sim()
    n = sim.mesh.n
    sim.driver.alpha = np.linspace(0.05, 0.9, n)
    assert _relative_error(sim) < 1e-5


def test_atomistic_jacobian_matches_finite_differences():
    """The atomistic driver defined sundials_jtn, which nothing ever called.

    Anisotropy is used rather than exchange because the atomistic exchange
    sets `jac = False`, so it is left out of the Jacobian: see
    test_interactions_marked_jac_false_are_left_out below.
    """
    mesh = CuboidMesh(nx=4, ny=4, nz=2)
    sim = AtomSim(mesh, name='jac_atom', integrator='cvode_bdf', use_jac=True)
    sim.set_mu_s(2 * mu_B)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 1.76e11
    sim.set_m(_random_unit_field(mesh.n, seed=3))
    sim.add(AtomAnisotropy(0.01 * k_B, axis=(0, 0, 1)))
    # the adaptive stabilisation term is not differentiated: see
    # test_adaptive_default_c_is_missing_from_the_jacobian
    sim.driver.default_c = 0
    assert _relative_error(sim) < 1e-5


def test_adaptive_default_c_is_missing_from_the_jacobian():
    """With default_c < 0 the right hand side and its Jacobian disagree.

    A negative default_c means the stabilisation term uses c = 6 |dm/dt|,
    which the right hand side applies and `llg_rhs_jtimes` then skips, its
    guard reading `if (default_c > 0)`. Differentiating it properly is not
    trivial, since c depends on m through dm/dt, so the term would have to be
    differentiated with c frozen at the value the right hand side used, which
    is available there as `fy`.

    It matters because -1 is the default of the atomistic driver, so every
    atomistic run with use_jac would carry it. This test records the gap.
    """
    mesh = CuboidMesh(nx=4, ny=4, nz=2)
    sim = AtomSim(mesh, name='jac_c', integrator='cvode_bdf', use_jac=True)
    sim.set_mu_s(2 * mu_B)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 1.76e11
    sim.set_m(_random_unit_field(mesh.n, seed=3))
    sim.add(AtomAnisotropy(0.01 * k_B, axis=(0, 0, 1)))

    sim.driver.default_c = 0
    assert _relative_error(sim) < 1e-5

    sim.driver.default_c = -1
    assert _relative_error(sim) > 1e-3, (
        'the adaptive stabilisation term now enters the Jacobian: if that '
        'was deliberate, this test should go')


def test_interactions_marked_jac_false_are_left_out():
    """Records what the `jac` flag currently does, which is a design question.

    `compute_effective_field_jac` sums only the interactions whose `jac` is
    True, so any other one is missing from the Jacobian. For a constant
    Zeeman field that is right, since it does not vary with m. For the
    demagnetising field it is not: the field is linear in m, so it does
    contribute to the derivative, and leaving it out makes the product
    approximate rather than exact. CVODE tolerates an approximate Jacobian,
    so this is defensible as a way of avoiding a transform per evaluation,
    but the flags are not consistent about it: the micromagnetic Demag sets
    False while the atomistic one sets True, and the micromagnetic exchange
    sets True while the atomistic one sets False.

    The number below is not a target. It is here so that changing a flag
    fails this test and brings someone back to the question.
    """
    from fidimag.micro import Demag as MicroDemag
    assert MicroDemag().jac is False
    assert AtomExchange(50 * k_B).jac is False

    sim = _micro_sim(with_demag=True)
    error = _relative_error(sim)
    assert error > 1e-3, (
        'the demagnetising field now enters the Jacobian: if that was '
        'deliberate, this test should go')


def test_use_jac_reaches_the_same_solution():
    """And the point of it all: it integrates, and agrees with the default."""
    def run(use_jac):
        mesh = CuboidMesh(3, 3, 10 / 3.0, 6, 6, 6, unit_length=1e-9)
        sim = MicroSim(mesh, name='jac_run_%d' % use_jac,
                       integrator='cvode_bdf', driver='llg', use_jac=use_jac)
        sim.Ms = 0.86e6
        sim.driver.alpha = 0.5
        sim.driver.set_tols(rtol=1e-10, atol=1e-10)
        sim.set_m((1, 0, 1))
        sim.add(UniformExchange(A=13e-12))
        sim.driver.run_until(5e-11)
        return sim.spin.copy()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.abs(run(True) - run(False)).max() < 1e-6
