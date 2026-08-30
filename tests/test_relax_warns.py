"""relax() must say so when the system has not relaxed (issue #118).

Reaching `max_steps` with `dmdt` still above `stopping_dmdt` used to return
quietly, and the result is indistinguishable from a relaxed state: the caller
gets a magnetisation either way. It now warns.
"""
import warnings

import numpy as np
import pytest

from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Zeeman

MU0 = 4 * np.pi * 1e-7


def _sim(name):
    mesh = CuboidMesh(nx=4, ny=4, nz=1, dx=2, dy=2, dz=2, unit_length=1e-9)
    sim = Sim(mesh, name=name)
    sim.set_Ms(8.0e5)
    sim.driver.alpha = 0.02
    sim.driver.gamma = 2.211e5
    sim.add(UniformExchange(A=1.3e-11))
    sim.add(Zeeman([0, 0, 0.05 / MU0]))
    sim.set_m((1, 0, 0))
    return sim


def test_warns_when_max_steps_is_reached():
    sim = _sim('relax_warn')
    with pytest.warns(RuntimeWarning, match='did not relax'):
        sim.driver.relax(max_steps=2, printing=False,
                         save_m_steps=None, save_vtk_steps=None)


def test_the_warning_reports_what_it_reached():
    """The message has to be actionable: both numbers, in the same units."""
    sim = _sim('relax_warn_msg')
    with pytest.warns(RuntimeWarning) as record:
        sim.driver.relax(max_steps=2, stopping_dmdt=0.01, printing=False,
                         save_m_steps=None, save_vtk_steps=None)
    message = str(record[0].message)
    assert 'max_steps=2' in message
    assert 'stopping_dmdt=0.01' in message
    assert 'max_dmdt=' in message


def test_silent_when_it_relaxes():
    sim = _sim('relax_quiet')
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        sim.driver.relax(max_steps=5000, stopping_dmdt=1.0, printing=False,
                         save_m_steps=None, save_vtk_steps=None)


def test_warns_when_the_loop_cannot_run():
    """`step` counts from the start of the simulation, not from this call.

    A second relax() with a max_steps below the current step does nothing at
    all, which is worth saying rather than returning as though it had worked.
    """
    sim = _sim('relax_noop')
    sim.driver.relax(max_steps=6, stopping_dmdt=1e-12, printing=False,
                     save_m_steps=None, save_vtk_steps=None)
    assert sim.driver.step >= 6
    with pytest.warns(RuntimeWarning, match='did nothing'):
        sim.driver.relax(max_steps=3, printing=False,
                         save_m_steps=None, save_vtk_steps=None)


def test_can_be_promoted_to_an_exception():
    sim = _sim('relax_raise')
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(RuntimeWarning):
            sim.driver.relax(max_steps=2, printing=False,
                             save_m_steps=None, save_vtk_steps=None)
