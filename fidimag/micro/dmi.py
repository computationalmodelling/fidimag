import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from energy import Energy
from fidimag.common.constant import mu_0
import fidimag.common.helper as helper


class DMI(Energy):

    """
        compute the DMI field in micromagnetics
    """

    def __init__(self, D, name='dmi', type='bulk'):
        """
        type could be 'interfacial' or 'bulk'
        """
        self.D = np.zeros(self.nxyz, dtype=np.float)
        self.D[:] = helper.init_scalar(D, self.mesh)
        self.name = name
        self.type = type

    def compute_field(self, t=0):
        if self.type == 'bulk':
            micro_clib.compute_dmi_field_bulk(self.spin,
                                              self.field,
                                              self.energy,
                                              self.Ms_inv,
                                              self.D,
                                              self.dx,
                                              self.dy,
                                              self.dz,
                                              self.nx,
                                              self.ny,
                                              self.nz,
                                              self.xperiodic,
                                              self.yperiodic)
        elif self.type == 'interfacial':
            micro_clib.compute_dmi_field_interfacial(self.spin,
                                                     self.field,
                                                     self.energy,
                                                     self.Ms_inv,
                                                     self.D,
                                                     self.dx,
                                                     self.dy,
                                                     self.dz,
                                                     self.nx,
                                                     self.ny,
                                                     self.nz,
                                                     self.xperiodic,
                                                     self.yperiodic)
        else:
            raise Exception(
                "Unsppourted dmi type:{}, avaiable type: 'bulk','interfacial'.".format(self.type))

        return self.field
