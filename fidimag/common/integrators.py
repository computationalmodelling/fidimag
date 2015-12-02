"""
Implement alternatives to sundials for testing and comparison purposes.

"""

def euler_step(t, y, h, f):
    """
    Numerical integration using the Euler method.

    Given the initial value problem y'(t) = f(t, y(t)), y(t_0) = y_0 one step
    of size h is y_{n+1} = y_n + h * f(t_n, y_n).

    """
    tp = t + h
    yp = y + h * f(t, y)
    return tp, yp


def runge_kutta_step(t, y, h, f):
    """
    Numerical integration using the classical Runge-Kutta method (RK4).

    Given the initial value problem y'(t) = f(t, y(t)), y(t_0) = y_0 one step
    of size h is y_{n+1} = y_n + h/6 * (k_1 + 2k_2 + 2k_3 + k4), where the
    weights are:
        k_1 = f(t_n,       y_n)
        k_2 = f(t_n + h/2, y_n + h/2 * k_1)
        k_3 = f(t_n + h/2, y_n + h/2 * k_2)
        k_4 = f(t_n + h,   y_n + h   * k_3).

    """
    k1 = f(t,           y)
    k2 = f(t + h / 2.0, y + h * k1 / 2.0)
    k3 = f(t + h / 2.0, y + h * k2 / 2.0)
    k4 = f(t + h,       y + h * k3)

    tp = t + h
    yp = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return tp, yp
