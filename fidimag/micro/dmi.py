import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from energy import Energy
from fidimag.common.constant import mu_0
import fidimag.common.helper as helper
import gc


class DMI(Energy):

    """
        compute the DMI field in micromagnetics
    """

    def __init__(self, D, name='dmi', type='bulk'):
        """
        type could be 'interfacial' or 'bulk'
        """
        self.D = D
        self.name = name
        self.jac = True
        self.type = type

    def setup(self, mesh, spin, Ms):
        super(DMI, self).setup(mesh, spin, Ms)
        self.Ds = np.zeros(self.n, dtype=np.float)
        self.Ds[:] = helper.init_scalar(self.D, self.mesh)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        if self.type == 'bulk':
            micro_clib.compute_dmi_field_bulk(m,
                                              self.field,
                                              self.energy,
                                              self.Ms_inv,
                                              self.Ds,
                                              self.dx,
                                              self.dy,
                                              self.dz,
                                              self.n,
                                              self.neighbours
                                              )

        elif self.type == 'interfacial':
            micro_clib.compute_dmi_field_interfacial(m,
                                                     self.field,
                                                     self.energy,
                                                     self.Ms_inv,
                                                     self.Ds,
                                                     self.dx,
                                                     self.dy,
                                                     self.dz,
                                                     self.n,
                                                     self.neighbours
                                                     )
        else:
            raise Exception(
                "Unsppourted dmi type:{}, avaiable type: 'bulk','interfacial'.".format(self.type))

        return self.field
