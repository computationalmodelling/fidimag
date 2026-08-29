"""The integrators exposed through SUNDIALS must agree with each other.

`cvode_bdf` is CVODE with BDF and a Newton iteration, `cvode_adams` is the
same solver with Adams and a fixed point iteration, and `arkode_dopri5` and
`arkode_rkf45` are explicit Runge-Kutta pairs from ARKODE. They are different methods, so
they take different steps, but integrating the same equation to the same
tolerance they have to land in the same place.
"""
import numpy as np
import pytest

import fidimag
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Zeeman

INTEGRATORS = ['cvode_bdf', 'cvode_adams', 'arkode_dopri5', 'arkode_rkf45',
               'arkode_dopri5_normalised', 'arkode_dopri5_openmp',
               'arkode_dopri5_normalised_openmp']


def _precession_sim(integrator):
    """A single macrospin precessing around an applied field."""
    mesh = CuboidMesh(nx=1, ny=1, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh, name='prec_%s' % integrator, integrator=integrator)
    sim.driver.set_tols(rtol=1e-10, atol=1e-10)
    sim.set_Ms(8.0e5)
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.211e5
    sim.add(Zeeman([0, 0, 0.1 / (4 * np.pi * 1e-7)]))
    sim.set_m((1, 0, 0))
    return sim


@pytest.mark.parametrize('integrator', INTEGRATORS)
def test_integrator_runs(integrator):
    sim = _precession_sim(integrator)
    sim.driver.run_until(1e-10)
    # the spin length is preserved, which the explicit methods do not do by
    # construction and only keep because the drivers normalise
    assert abs(np.linalg.norm(sim.spin) - 1.0) < 1e-8
    steps, nfevals, _ = sim.driver.stat()
    assert steps > 0 and nfevals > 0


def test_integrators_agree():
    """Same trajectory, to a tolerance well above the requested 1e-10."""
    ts = np.linspace(0, 2e-10, 11)[1:]
    traj = {}
    for integrator in INTEGRATORS:
        sim = _precession_sim(integrator)
        out = []
        for t in ts:
            sim.driver.run_until(t)
            out.append(sim.spin.copy())
        traj[integrator] = np.array(out)

    ref = traj['cvode_bdf']
    for integrator in INTEGRATORS[1:]:
        deviation = np.abs(traj[integrator] - ref).max()
        assert deviation < 1e-6, \
            '{} deviates from cvode_bdf by {:.2e}'.format(integrator, deviation)


def test_normalised_variant_conserves_the_spin_length_better():
    """`_normalised` rescales the spins after every accepted step.

    The length of m is otherwise kept near one by the `c * (1 - m^2) * m`
    correction term of the LLG right hand side, which makes |m| = 1 an
    attracting solution rather than imposing it. Projecting as well is more
    accurate, and should cost close to nothing: the projection is O(n)
    arithmetic once per step, against the effective field evaluations, with
    their demagnetising transforms, that each step already pays for.
    """
    plain = _precession_sim('arkode_dopri5')
    plain.driver.run_until(2e-10)

    projected = _precession_sim('arkode_dopri5_normalised')
    projected.driver.run_until(2e-10)

    err_plain = abs(np.linalg.norm(plain.spin) - 1.0)
    err_projected = abs(np.linalg.norm(projected.spin) - 1.0)
    assert err_projected < err_plain

    # and it is the same trajectory
    assert np.abs(projected.spin - plain.spin).max() < 1e-6

    # for the same work, to within the step controller's discretion
    steps_plain = plain.driver.stat()[0]
    steps_projected = projected.driver.stat()[0]
    assert steps_projected <= 1.2 * steps_plain


def test_openmp_matches_serial():
    """The OpenMP N_Vector only threads the integrator's own arithmetic.

    It changes how the vector operations are executed, not what they compute,
    so the two must follow the same trajectory step for step.
    """
    serial = _precession_sim('arkode_dopri5')
    serial.driver.run_until(2e-10)

    threaded = _precession_sim('arkode_dopri5_openmp')
    threaded.driver.run_until(2e-10)

    assert np.abs(threaded.spin - serial.spin).max() < 1e-9
    assert threaded.driver.stat()[0] == serial.driver.stat()[0]


def test_unknown_butcher_table_is_rejected():
    from fidimag.common.integrators import ErkSolver

    def rhs(t, y, ydot):
        ydot[:] = -y[:]
        return 0

    with pytest.raises(RuntimeError):
        ErkSolver(np.ones(3), rhs, table='ARKODE_NOT_A_REAL_TABLE')


def test_unknown_integrator_is_rejected():
    mesh = CuboidMesh(nx=1, ny=1, nz=1)
    with pytest.raises(NotImplementedError):
        Sim(mesh, name='bad', integrator='not_an_integrator')


@pytest.mark.parametrize('legacy,canonical', [
    ('sundials', 'cvode_bdf'),
    ('sundials_diag', 'cvode_bdf_diag'),
    ('sundials_openmp', 'cvode_bdf_openmp'),
    ('sundials_diag_openmp', 'cvode_bdf_diag_openmp'),
])
def test_pre_4_0_names_still_work_with_a_warning(legacy, canonical):
    """Scripts written before the rename must keep running."""
    from fidimag.common.driver_base import DriverBase

    with pytest.warns(DeprecationWarning, match=canonical):
        assert DriverBase._canonical_integrator(legacy) == canonical

    # and the canonical name itself warns about nothing
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert DriverBase._canonical_integrator(canonical) == canonical


def test_legacy_name_runs_a_simulation():
    mesh = CuboidMesh(nx=1, ny=1, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    with pytest.warns(DeprecationWarning):
        sim = Sim(mesh, name='legacy', integrator='sundials')
    sim.driver.set_tols(rtol=1e-10, atol=1e-10)
    sim.set_Ms(8.0e5)
    sim.add(Zeeman([0, 0, 0.1 / (4 * np.pi * 1e-7)]))
    sim.set_m((1, 0, 0))
    sim.driver.run_until(1e-11)
    assert abs(np.linalg.norm(sim.spin) - 1.0) < 1e-8
