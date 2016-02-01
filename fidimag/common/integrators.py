"""
Implement the time integrators used in the simulation of magnetisation dynamics.

"""
import fidimag.extensions.cvode as cvode

EPSILON = 1e-14


class SundialsIntegrator(object):
    def __init__(self, spins, rhs):
        self.cvode = cvode.CvodeSolver(spins, rhs)
        self.set_tols()

    @property
    def rhs_evals(self):
        """
        This function tries to get the values from the CVODE statistics. For a
        'sim' simulation object, this is done starting from calling
        sim.vode.stats()

        According to the CVODE version, this call can generate a string:

        CvodeSolver(nsteps = 18,
                    nfevals = 32,
                    njevals = 14.
                    )

        where:

        nsteps  --> number of steps taken by CVODE
        nfevals --> number of calls to the user's f function
                    (I guess this is what we need)
        njevals --> the cumulative number of calls to the Jacobian function

        So, for example,  we can regex search any number preceded by
            "nfevals = "
        to get the number of evaluations of the RHS and convert the
        result to an integer

        OR it can give a tuple with 3 values, which must be in the same
        order than before

        For now, we are only interested in the RHS evaluations, so we
        return a single value

        """
        stat = self.stat()

        if isinstance(stat, str):
            out = int(re.search(r'(?<=nfevals\s=\s)[0-9]*', stat).group(0))
        elif isinstance(stat, tuple):
            out = stat[1]
        else:
            raise NotImplementedError('Cannot retrieve the values'
                                      'from CVODE stats')

        return out

    def run_until(self, t):
        return self.cvode.run_until(t)

    def reset(self, spins, t):
        self.cvode.reset(spins, t)

    def set_tols(self, rtol=1e-8, atol=1e-10):
        self.cvode.set_options(rtol, atol)

    def set_initial_value(self, spins, t, reuse_memory=1):
        self.cvode.set_initial_value(spins, t, reuse_memory)

    def stat(self):
        return self.cvode.stat()

    def get_current_step(self):
        return self.cvode.get_current_step()

    @property
    def y(self):
        return self.cvode.y


class StepIntegrator(object):
    def __init__(self, spins, rhs, step="euler"):
        self.spins = spins
        self.rhs = rhs
        self.t = 0
        self.h = 1e-14
        self.steps = 0
        if step == "euler":
            self.single_step = euler_step
        elif step == "rk4":
            self.single_step = runge_kutta_step
        else:
            raise NotImplemented("step must be euler or rk4")

    @property
    def rhs_evals(self):
        return self.steps

    def run_until(self, t):
        while abs(self.t - t) > EPSILON:
            self.t, self.spins = self.single_step(self.t, self.spins, self.h, self.rhs)
            self.steps += 1
            if self.t > t:
                break
        return 0

    def reset(self, spins, t):
        self.spins = spins
        self.t = t
        self.steps = 0

    # same methods as SundialsIntegrator below
    def set_tols(self, rtol=1e-8, atol=1e-10):  
        pass

    def set_initial_value(self, spins, t, reuse_memory=1):
        self.spins = spins
        self.t = t

    def stat(self):
        return self.steps

    def get_current_step(self):
        return self.h

    @property
    def y(self):
        return self.spins


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
