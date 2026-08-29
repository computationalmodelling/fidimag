"""muMAG standard problem 4, against the OOMMF reference solution.

The reference is `examples/micromagnetic/std4/oommf.txt`, an OOMMF run of the
same problem sampled every picosecond out to a nanosecond. The comparison
lived only in that example, so nothing ran it; this is the same check as a
test.

The problem is the 500 x 125 x 3 nm permalloy film of the muMAG standard
problem 4, relaxed to its S-state and then reversed by field 1,
mu0 * H = (-24.6, 4.3, 0) mT.
"""
import os

import numpy as np
import pytest

from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Demag, Zeeman

MU0 = 4 * np.pi * 1e-7
A = 1.3e-11
MS = 8.0e5
GAMMA = 2.211e5

REFERENCE = os.path.join(os.path.dirname(__file__), os.pardir, 'examples',
                         'micromagnetic', 'std4', 'oommf.txt')


def _mesh():
    return CuboidMesh(nx=200, ny=50, nz=1, dx=2.5, dy=2.5, dz=3,
                      unit_length=1e-9)


def _s_state(mesh):
    """Relax to the S-state that the problem starts from."""
    sim = Sim(mesh, name='std4_relax')
    sim.driver.set_tols(rtol=1e-10, atol=1e-10)
    sim.driver.alpha = 0.5
    sim.driver.gamma = GAMMA
    sim.Ms = MS
    sim.driver.do_precession = False
    sim.set_m((1, 0.25, 0.1))
    sim.add(UniformExchange(A=A))
    sim.add(Demag())
    sim.driver.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
                     save_m_steps=None, save_vtk_steps=None)
    return sim.spin.copy()


def _reversal(mesh, m0, ts, integrator='cvode_bdf'):
    sim = Sim(mesh, name='std4_dyn', integrator=integrator)
    sim.driver.set_tols(rtol=1e-10, atol=1e-10)
    sim.driver.alpha = 0.02
    sim.driver.gamma = GAMMA
    sim.Ms = MS
    sim.set_m(m0)
    sim.add(UniformExchange(A=A))
    sim.add(Demag())
    sim.add(Zeeman([-24.6e-3 / MU0, 4.3e-3 / MU0, 0], name='H'))

    out = []
    for t in ts:
        sim.driver.run_until(t)
        out.append(sim.spin.reshape(-1, 3).mean(axis=0).copy())
    return np.array(out)


@pytest.mark.slow
def test_std_problem_4_matches_oommf():
    reference = np.loadtxt(REFERENCE)
    # every tenth sample: the deviation is measured over the whole nanosecond
    # either way, and the integration, not the sampling, is what costs
    ts, m_oommf = reference[::10, 0], reference[::10, 1:4]

    mesh = _mesh()
    m = _reversal(mesh, _s_state(mesh), ts)

    deviation = np.abs(m - m_oommf)
    print('deviation from OOMMF: max {:.3e}, rms {:.3e}, last {}'.format(
        deviation.max(), np.sqrt((deviation ** 2).mean()), deviation[-1]))

    # Measured 3.53e-05 and 1.31e-05 when this was written. The thresholds
    # leave room for the relaxed S-state to shift slightly, which moves the
    # whole trajectory, without admitting a real regression
    assert deviation.max() < 1e-4
    assert np.sqrt((deviation ** 2).mean()) < 5e-5

    # the reversal happens: <m_x> starts positive and ends negative
    assert m[0, 0] > 0.9
    assert m[-1, 0] < -0.9


@pytest.mark.slow
def test_std_problem_4_agrees_between_integrators():
    """The explicit method must reproduce the reference solution too."""
    reference = np.loadtxt(REFERENCE)
    ts, m_oommf = reference[::50, 0], reference[::50, 1:4]

    mesh = _mesh()
    m0 = _s_state(mesh)

    implicit = _reversal(mesh, m0, ts, integrator='cvode_bdf')
    explicit = _reversal(mesh, m0, ts, integrator='arkode_dopri5')

    assert np.abs(explicit - implicit).max() < 1e-6
    assert np.abs(explicit - m_oommf).max() < 1e-4
