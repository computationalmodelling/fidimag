import numpy as np
import fidimag.extensions.baryakhtar_clib as clib
from .micro_driver import MicroDriver
from .relax import Relaxation
from .relax import Laplace


class LLBarFull(MicroDriver):

    def __init__(self, mesh, spin, Ms, field, alpha, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False,
                 chi=1e-3
                 ):

        # Inherit from the driver class
        super(LLBarFull, self).__init__(mesh, spin, Ms, field,
                                        alpha, pins, interactions, name,
                                        data_saver,
                                        integrator=integrator,
                                        use_jac=False
                                        )

        self.chi = chi
        self.lap = Laplace(mesh)
        self.add(Relaxation(chi))

        self.beta = 0

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)
        delta_h = self.lap.compute_laplace_field(self.field, self._Ms)

        clib.compute_llg_rhs_baryakhtar(ydot,
                                        self.spin,
                                        self.field,
                                        delta_h,
                                        self.alpha,
                                        self.beta,
                                        self._pins,
                                        self.gamma,
                                        self.n,
                                        self.do_precession)

        #ydot[:] = self.dm_dt[:]

        return 0


class LLBar(MicroDriver):

    def __init__(self, mesh, spin, Ms, field, alpha, pins,
                 interactions,
                 name,
                 data_saver,
                 integrator='sundials',
                 use_jac=False,
                 ):

        # Inherit from the driver class
        super(LLBar, self).__init__(mesh, spin, Ms, field,
                                    alpha, pins, interactions, name,
                                    data_saver,
                                    integrator='sundials',
                                    use_jac=False
                                    )

        self.lap = Laplace(mesh)

        self.field_perp = np.zeros(3 * self.n, dtype=np.float)

        self.beta = 0

    def sundials_rhs(self, t, y, ydot):

        self.t = t

        # already synchronized when call this funciton
        # self.spin[:]=y[:]

        self.compute_effective_field(t)
        clib.compute_perp_field(
            self.spin, self.field, self.field_perp, self.n)
        delta_h = self.lap.compute_laplace_field(self.field_perp, self._Ms)

        clib.compute_llg_rhs_baryakhtar_reduced(ydot,
                                                self.spin,
                                                self.field,
                                                delta_h,
                                                self.alpha,
                                                self.beta,
                                                self._pins,
                                                self.gamma,
                                                self.n,
                                                self.do_precession,
                                                self.default_c)

        #ydot[:] = self.dm_dt[:]

        return 0
