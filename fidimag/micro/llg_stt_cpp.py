from __future__ import division

import types

import fidimag.extensions.clib as clib
import numpy as np

from .llg import LLG
import fidimag.common.helper as helper


class LLG_STT_CPP(LLG):

    def __init__(self, mesh, name='unnamed'):
        """Simulation object.

        *Arguments*

          name : the Simulation name (used for writing data files, for examples)

        """
        super(LLG_STT_CPP, self).__init__(mesh, name=name)

        self._p = np.zeros(3 * self.n, dtype=np.float)
        self._a_J = np.zeros(self.n, dtype=np.float)

        self.a_J = 1
        self.beta = 0
        self.J_time_fun = None
        self.p = (0,0,1)

    def get_p(self):
        return self._p

    def set_p(self, value):
        self._p[:] = helper.init_vector(value, self.mesh)
    p = property(get_p, set_p)

    def get_a_J(self):
        return self._a_J

    def set_a_J(self, value):
        self._a_J[:] = helper.init_scalar(value, self.mesh)

    a_J = property(get_a_J, set_a_J)

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        if self.J_time_fun is not None:
            clib.compute_llg_stt_cpp(ydot,
                                 self.spin,
                                 self.field,
                                 self._p,
                                 self.alpha,
                                 self._pins,
                                 self._a_J*self.J_time_fun(t),
                                 self.beta,
                                 self.gamma,
                                 self.n)
        else:
            clib.compute_llg_stt_cpp(ydot,
                                 self.spin,
                                 self.field,
                                 self._p,
                                 self.alpha,
                                 self._pins,
                                 self._a_J,
                                 self.beta,
                                 self.gamma,
                                 self.n)
