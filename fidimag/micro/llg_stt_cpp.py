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

        self.a_J = 1
        self.beta = 0

    def get_p(self):
        return self._p

    def set_p(self, value):
        self._p[:] = helper.init_vector(value, self.mesh)

    p = property(get_p, set_p)

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)

        if isinstance(self.a_J,types.FunctionType):
            _a_J = self.a_J(t)
        else:
            _a_J = self.a_J
        
        clib.compute_llg_stt_cpp(ydot,
                                 self.spin,
                                 self.field,
                                 self._p,
                                 self.alpha,
                                 self._pins,
                                 self.beta,
                                 _a_J,
                                 self.gamma,
                                 self.n)
