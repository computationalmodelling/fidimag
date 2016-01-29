"""
Sanity testing of the integrators in integrators.py.

"""
import matplotlib as mpl
mpl.use('Agg')
import pytest
import numpy as np
from math import ceil
from integrators import euler_step, runge_kutta_step, StepIntegrator


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
        tp, yp = integrator(t, ys[i], stepsize, f)
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
    print ts
    print ts[1:]
    ys = np.zeros(ts.shape[0])

    ys[0] = y_true(ts[0])  # known initial value

    integrator = StepIntegrator(ys[0], f, step="euler")
    integrator.h = stepsize_internal
    for i, t in enumerate(ts[1:]):
        assert integrator.run_until(t) == 0
        ts[i+1] = integrator.t
        ys[i+1] = integrator.y
    assert 85 < ys[-1] < 100
    return ts, ys

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
