import numpy as np


class DriverBase(object):

    """
    The DriverBase class is a supper class from which specific drivers can be implemented.

    In Fidimag, a driver implements the mechanism to drive the magnetization to a new state.
    The behind mechanisms could be LLG-based equations, Energy minimization, Monte Carlo or 
    other interesting rules.

    Each derived class should override the next_step method.

    """

    def __init__(self, mesh, spin, magnitude, pins, interactions, field):
        self.t = 0
        self.step = 0
        self.mesh = mesh
        self.spin = spin
        self.n = mesh.n
        #magnitude could be Ms or mu_s.
        self._magnitude = magnitude 
        self._pins = pins
        self.interactions = interactions
        self.field = field
        self.previous_spin = np.copy(spin)
    
    def next_step(self, time=None):

        raise Exception('Derived classes should override this method!')
