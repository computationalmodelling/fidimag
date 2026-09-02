import fidimag
import fidimag.common.constant as C
import numpy as np


Ms = 0.86e6
A = 13e-12
Ku = 0.4e6

# Scaled energy of the relaxed 1D domain wall. Note that the wall is free to
# slide along x at no energy cost, so at an aggressive `tmax` it can end up
# translated: the energy identifies the minimum, the profile position does not
E_DW_REFERENCE = 1.9603112400e-01


def _setup_1D_DW(tmax=None, driver='steepest_descent'):
    """1D micromagnetic sample with a domain wall as the initial state"""

    mesh = fidimag.common.CuboidMesh(nx=100, ny=1, nz=1, dx=1, dy=1, dz=1,
                                     unit_length=1e-9)

    sim = fidimag.micro.Sim(mesh, name='1Dmicro_sd', driver=driver)

    sim.set_Ms(Ms)
    sim.add(fidimag.micro.UniformExchange(A))
    sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))

    xs = mesh.coordinates[:, 0]
    centre_x = (xs.max() + xs.min()) * 0.5 + xs.min()

    def m_initial(r):
        if r[0] < centre_x:
            return (0, 0.1, -.9)
        else:
            return (0, 0.1, .9)
    sim.set_m(m_initial)

    mesh_vol = mesh.n * mesh.dx * mesh.dy * mesh.dz * 1e-27
    sim.driver.energyScale = C.mu_0 * (Ms ** 2) * 0.5 * mesh_vol
    if tmax is not None:
        sim.driver.tmax = tmax

    return sim


def _total_energy(sim):
    return (sum(obj.compute_energy() for obj in sim.driver.interactions)
            / sim.driver.energyScale)


def _dw_MAE(sim):
    """Mean absolute error against the analytical wall, within the DW width"""

    mz = sim.spin.reshape(-1, 3)[:, 2]
    x = sim.mesh.coordinates[:, 0]
    y = np.tanh((x * 1e-9 - 50 * 1e-9) / np.sqrt(A / Ku))

    deltaB = np.sqrt(A / Ku) * 1e9
    xc = 50
    ftr = np.logical_and(x <= xc + deltaB, x >= xc - deltaB)

    return np.mean(np.abs(y[ftr] - mz[ftr]))


def test_steepest_descent_1D_DW():

    sim = _setup_1D_DW()
    sim.driver.minimise(stopping_dm=1e-9, max_steps=3000, printing=False)

    MAE_dw = _dw_MAE(sim)
    print('Domain wall MAE', MAE_dw)
    assert MAE_dw < 0.02


def test_steepest_descent_tau_is_updated():
    """The Barzilai-Borwein step size must reach Python from the C routine

    `sd_compute_step` takes `tau` by value, so a step size computed there and
    assigned to the parameter is lost: it has to be returned and assigned back
    on the Python side. When that is missed, `self.tau` keeps the value of
    `initial_t_step` forever and the whole BB scheme silently degrades into
    fixed-step descent.
    """

    initial_t_step = 1e-2
    sim = _setup_1D_DW()
    # stopping_dm=0 so that the run cannot stop before max_steps
    sim.driver.minimise(stopping_dm=0.0, max_steps=20, printing=False,
                        initial_t_step=initial_t_step)

    assert sim.driver.tau != initial_t_step


def test_steepest_descent_tau_stays_positive():
    """A non-positive BB quotient must not be turned into an uphill step

    The quotient is negative where the curvature along the last step is
    non-positive; taking it with its sign steps against -m x m x H.
    """

    sim = _setup_1D_DW()
    taus = []
    run_step = sim.driver.run_step_CLIB

    def record():
        run_step()
        taus.append(sim.driver.tau)
    sim.driver.run_step_CLIB = record

    sim.driver.minimise(stopping_dm=0.0, max_steps=200, printing=False)

    taus = np.array(taus)
    assert (taus > 0.0).all()
    assert (taus <= sim.driver.tmax).all()


def test_steepest_descent_energy_guard_large_tmax():
    """The energy guard is what makes an aggressive step ceiling usable

    With a large `tmax` the BB quotient can propose a step that overshoots
    badly. Unguarded, that step is taken regardless: the iteration then either
    spends most of its steps recovering or runs away altogether, and because
    the BB quotient is accumulated with an OpenMP reduction (whose summation
    order is not fixed) which of the two happens is not even reproducible from
    run to run. The non-monotone sufficient-decrease test backtracks those
    steps instead, which both bounds the damage and makes the trajectory
    repeatable.

    Only the guarded run is asserted on here: the unguarded one fails
    intermittently by construction, so testing *its* failure would be flaky.
    """

    sim = _setup_1D_DW(tmax=10.0)
    sim.driver.minimise(stopping_dm=1e-9, max_steps=3000, printing=False,
                        initial_t_step=1e-2, energy_guard=True)

    # Converged to the domain wall, and not to some overshot configuration
    assert sim.driver.step < 3000
    assert abs(_total_energy(sim) - E_DW_REFERENCE) < 1e-6
    # The guard is doing something here: trial steps really were rejected
    assert sim.driver.nEval > sim.driver.step + 1


def test_steepest_descent_energy_guard_matches_default():
    """At the default (conservative) tmax the guard must not change the answer

    There the BB steps rarely trip the acceptance test, so the guard should be
    a no-op on the result and only cost the energy evaluations.
    """

    sim_plain = _setup_1D_DW()
    sim_plain.driver.minimise(stopping_dm=1e-9, max_steps=3000, printing=False)

    sim_guard = _setup_1D_DW()
    sim_guard.driver.minimise(stopping_dm=1e-9, max_steps=3000, printing=False,
                              energy_guard=True)

    assert _dw_MAE(sim_guard) < 0.02
    assert abs(_total_energy(sim_guard) - _total_energy(sim_plain)) < 1e-6


if __name__ == "__main__":
    test_steepest_descent_1D_DW()
    test_steepest_descent_tau_is_updated()
    test_steepest_descent_tau_stays_positive()
    test_steepest_descent_energy_guard_large_tmax()
    test_steepest_descent_energy_guard_matches_default()


def test_minimisers_are_quiet_unless_asked(capsys):
    """
    No minimiser prints: progress goes to the `fidimag` logger.

    The per step lines are at DEBUG, the reason the iteration stopped at INFO
    and a failure to converge at WARNING, so a run that converges is silent
    unless that logger is turned up.
    """
    for driver, kwargs in [('steepest_descent', dict(stopping_dm=1e-4)),
                           ('hubert_minimiser', dict(stepControl='BB',
                                                     mXgradE_tol=1e-2,
                                                     stopping_dE=-1.0))]:
        sim = _setup_1D_DW(driver=driver)
        sim.driver.minimise(max_steps=200, **kwargs)
        out = capsys.readouterr()
        assert out.out == '', '%s wrote to stdout: %r' % (driver, out.out[:120])
        assert out.err == '', '%s wrote to stderr: %r' % (driver, out.err[:120])


def test_stopping_torque_stops_on_the_residual():
    """
    `stopping_torque` tests the residual, `stopping_dm` a displacement.

    `max_dm` is the product of the step length and the torque, and
    Barzilai-Borwein step lengths swing over orders of magnitude, so it can
    fall below its threshold on a short step while the residual is unchanged.
    """
    sim = _setup_1D_DW()
    sim.driver.minimise(stopping_torque=1e-4, max_steps=20000, printing=False)

    assert sim.driver.max_torque() < 1e-4
    assert _dw_MAE(sim) < 0.02

    # the displacement criterion, asked for the same tightness, stops with a
    # residual orders of magnitude larger
    loose = _setup_1D_DW()
    loose.driver.minimise(stopping_dm=1e-4, max_steps=20000, printing=False)
    assert loose.driver.max_torque() > 1e-4


def test_energy_guard_is_what_allows_a_large_tmax():
    """
    Above a certain step ceiling the guard is the only thing keeping the
    iteration on the right configuration.

    The collapsed wall is stationary, its torque satisfying the same
    criterion, so no stopping test can catch it: the energy has to be looked
    at while the step is taken.
    """
    unguarded = _setup_1D_DW(tmax=10.0)
    unguarded.driver.minimise(stopping_torque=1e-4, max_steps=20000,
                              printing=False, energy_guard=False)
    assert unguarded.driver.max_torque() < 1e-4, 'it did converge, but ...'
    assert _dw_MAE(unguarded) > 0.5, '... to the wrong configuration'

    guarded = _setup_1D_DW(tmax=10.0)
    guarded.driver.minimise(stopping_torque=1e-4, max_steps=20000,
                            printing=False, energy_guard=True)
    assert guarded.driver.max_torque() < 1e-4
    assert _dw_MAE(guarded) < 0.02
