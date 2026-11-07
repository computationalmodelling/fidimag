"""The fixed step integrators of the drivers: euler, rk4 and their _openmp names."""
import warnings

import numpy as np
import pytest

import fidimag.common.constant as const
from fidimag.common import CuboidMesh
from fidimag.common.integrators import runge_kutta_step


def micro_sim(integrator, name='steps'):
    from fidimag.micro import Sim, UniformExchange, Zeeman
    sim = Sim(CuboidMesh(nx=4, ny=1, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9), name=name, integrator=integrator)
    sim.set_Ms(8e5)
    sim.alpha, sim.gamma = 0.5, 2.211e5
    sim.set_m(lambda r: (np.sin(3e8 * r[0]), np.cos(2e8 * r[0]), 0.3))
    sim.add(UniformExchange(1.3e-11))
    sim.add(Zeeman((1e4, 0, 3e4)))
    return sim


def run(integrator, stepsize=None, t=2e-12, sim=None):
    sim = sim or micro_sim(integrator)
    if stepsize is not None:
        sim.driver.integrator.stepsize = stepsize
    sim.driver.run_until(t)
    return sim.spin.copy()


def test_rk4_stages_do_not_share_an_array():
    """A right hand side that returns the same array every time (as `step_rhs` of the drivers does) must not make the four stages one."""
    buffer = np.zeros(2)

    def f(t, y):
        buffer[:] = (-y[0], y[0] - 2 * y[1] + t)
        return buffer

    def g(t, y):
        return np.array([-y[0], y[0] - 2 * y[1] + t])

    y0 = np.array([1.0, 0.5])
    shared = runge_kutta_step(0.1, y0.copy(), 0.05, f)
    fresh = runge_kutta_step(0.1, y0.copy(), 0.05, g)
    assert shared[0] == fresh[0] and np.array_equal(shared[1], fresh[1])


def test_euler_is_first_order():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        reference = run('cvode_bdf')
        coarse = np.abs(run('euler', 2e-14) - reference).max()
        fine = np.abs(run('euler', 1e-14) - reference).max()
    assert coarse / fine == pytest.approx(2, rel=0.2)


def test_rk4_follows_the_reference():
    """(Its four stages were one, so it was no better than Euler: an error of 1e-2 here, against 1e-6 or less)"""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        reference = run('cvode_bdf')
        assert np.abs(run('rk4', 1e-14) - reference).max() < 1e-5


@pytest.mark.parametrize('name', ['euler', 'rk4'])
def test_the_openmp_names_give_the_same_steps(name):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.array_equal(run(name + '_openmp', 1e-14), run(name, 1e-14))


@pytest.mark.parametrize('integrator', ['euler', 'rk4'])
def test_set_m_restarts_the_integrator_from_the_new_state(integrator):
    """`Sim.set_m` tells the integrator the new state with `set_initial_value`, which did nothing for a step integrator, so that after a run it
    went on from where it had been."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sim = micro_sim(integrator)
        sim.driver.integrator.stepsize = 1e-14
        sim.driver.run_until(1e-12)
        new = lambda r: (np.cos(2e8 * r[0]), 0.4, np.sin(3e8 * r[0]))
        sim.set_m(new)
        t = sim.driver.t
        sim.driver.run_until(t + 1e-12)
        restarted = sim.spin.copy()

        fresh = micro_sim(integrator, name='steps_fresh')
        fresh.set_m(new)
        fresh.driver.integrator.stepsize = 1e-14
        fresh.driver.run_until(1e-12)
    assert np.abs(restarted - fresh.spin).max() < 1e-12


@pytest.mark.parametrize('integrator', ['euler', 'rk4'])
def test_atomistic_drivers_have_the_step_integrators(integrator):
    from fidimag.atomistic import Sim, UniformExchange, Zeeman

    def atomistic(name):
        sim = Sim(CuboidMesh(nx=4, ny=1, nz=1), name='atomistic_' + name, integrator=name)
        sim.driver.set_default_options(gamma=const.gamma)
        sim.mu_s = const.mu_B
        sim.alpha = 0.5
        sim.set_m(lambda r: (np.sin(3 * r[0]), np.cos(2 * r[0]), 0.3))
        sim.add(UniformExchange(1e-21))
        sim.add(Zeeman((0, 0, 1.0)))
        return sim

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sim, reference = atomistic(integrator), atomistic('cvode_bdf')
        sim.driver.integrator.stepsize = 1e-16
        sim.driver.run_until(5e-14)
        reference.driver.run_until(5e-14)
    assert np.abs(sim.spin - reference.spin).max() < 5e-3
