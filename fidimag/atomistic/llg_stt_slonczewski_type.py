from __future__ import division

import types

import fidimag.extensions.clib as clib
import numpy as np

from .llg import LLG
import fidimag.common.helper as helper


class LLG_STT_Slonczewski(LLG):

    def __init__(self, mesh, name='unnamed'):
        """Simulation object.

        *Arguments*

          name : the Simulation name (used for writing data files, for examples)

        """
        super(LLG_STT_Slonczewski, self).__init__(mesh, name=name)

        self._p = np.zeros(3 * self.n, dtype=np.float)

        self.u0 = 1

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

        if isinstance(self.u0,types.FunctionType):
            _u = self.u0(t)
        else:
            _u = self.u0
        
        clib.compute_llg_stt_slonczewski_type(ydot,
                                 self.spin,
                                 self.field,
                                 self._p,
                                 self.alpha,
                                 _u,
                                 self.gamma,
                                 self.n)
