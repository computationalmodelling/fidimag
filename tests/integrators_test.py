"""
Sanity testing of the integrators in integrators.py.

"""
from __future__ import print_function
import matplotlib as mpl
mpl.use('Agg')
import pytest
import numpy as np
from math import ceil
from .integrators import euler_step, runge_kutta_step, StepIntegrator, ScipyIntegrator


interval = (0, 10)
y_true = lambda t: t ** 2
f = lambda t, y: 2 * t


@pytest.mark.parametrize("integrator,stepsize", [
    (euler_step, 0.2),
    (euler_step, 0.5),
    (runge_kutta_step, 0.2),
    (runge_kutta_step, 0.5)])
def test_step(integrator, stepsize, debug=False):
    ts = np.arange(interval[0], interval[1], step=stepsize)
    ys = np.zeros(ts.shape[0])

    ys[0] = y_true(ts[0])  # known initial value

    for i, t in enumerate(ts[:-1]):
        tp, yp, evals = integrator(t, ys[i], stepsize, f)
        ts[i+1] = tp
        ys[i+1] = yp
    if not debug:
        assert 85 < ys[-1] < 100
    return ts, ys

@pytest.mark.parametrize("integrator,stepsize_reported,stepsize_internal", [
    ("euler", 0.2, 0.2),
    ("euler", 0.2, 0.05),
    ("rk4", 0.2, 0.2),
    ("rk4", 0.2, 0.05)])
def test_step_integrator(integrator, stepsize_reported, stepsize_internal):
    ts = np.arange(interval[0], interval[1], step=stepsize_reported)
    print(ts)
    print(ts[1:])
    ys = np.zeros(ts.shape[0])

    ys[0] = y_true(ts[0])  # known initial value

    integrator = StepIntegrator(ys[0], f, step="euler", stepsize=stepsize_internal)
    for i, t in enumerate(ts[1:]):
        integrator.run_until(t)
        ts[i+1] = integrator.t
        ys[i+1] = integrator.y
    assert 85 < ys[-1] < 100
    return ts, ys


def test_scipy_integrator():
    y_true = lambda t: np.sin(t) + t
    f = lambda t, y: np.cos(t) + 1  # derivative of f we'll use for integration

    ts = np.arange(interval[0], interval[1], step=1.0)
    ys = np.zeros(ts.shape[0])
    ys[0] = y_true(ts[0])  # known initial value

    y = np.zeros(1)  # scipy wants an np.array
    y[0] = ys[0]
    od = ScipyIntegrator(y, f)
    od.set_options(rtol=1e-5, atol=1e-5)

    for i, t in enumerate(ts[1:]):
        od.run_until(t)
        ts[i+1] = od.t
        ys[i+1] = od.y
        print("t = {}".format(t))


    print(len(od.internal_timesteps))
    print(od.rhs_evals)
    return y_true, ts, ys, od.internal_timesteps


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    ts_fine = np.linspace(interval[0], interval[1], num=100)
    plt.plot(ts_fine, y_true(ts_fine), label="y=x**2")

    for h in (0.5, 1):
        plt.plot(*test_step(euler_step, h, debug=True), marker="o", linestyle="dashed", label="euler h={}".format(h))
        plt.plot(*test_step(runge_kutta_step, h, debug=True), marker="o", linestyle="dashed", label="RK4 h={}".format(h))
    plt.plot(*test_step_integrator("euler", 0.2, 0.05), marker="o", linestyle="dashed", label="euler h=0.05 int")
    plt.legend(loc=0)
    plt.savefig("test_integrators.png")
    plt.clf()

    y_true, ts, ys, ts_internal = test_scipy_integrator()
    plt.plot(ts_fine, y_true(ts_fine), label="y=sin(x)+x")
    plt.plot(ts, ys, "go", label="dopri reported")
    for t in ts_internal[:-1]:
        plt.plot((t, t), (0, 1), 'r-')
    plt.plot((ts_internal[-1], ts_internal[-1]), (0, 1), 'r-', label="internal timesteps")
    plt.legend(loc=0)
    plt.savefig("test_integrators_dopri.png")
