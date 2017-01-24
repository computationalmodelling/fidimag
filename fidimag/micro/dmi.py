import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from .energy import Energy
from fidimag.common.constant import mu_0
import fidimag.common.helper as helper
import gc


class DMI(Energy):

    """

    Compute the Dzyaloshinskii-Moriya interaction in the micromagnetic
    framework

    ARGUMENTS: ----------------------------------------------------------------

    D       ::

               DMI vector norm which can be specified as an int, float, (X * n)
               array (X=6 for bulk DMI and X=4 for interfacial DMI), (n) array
               or spatially dependent scalar field function.

               int, float: D will have the same magnitude for every NN of the
               spins at every mesh node, given by this magnitude

               (n) array: D will have the same magnitude for every NN on every
               mesh node but will vary according to the values of the array,
               i.e. the value of D for the 6 NNs at the i-th mesh site will be
               given by the i-th value of the array

               (X * n) array: Manually specify the DMI vector norm for every NN
               at every mesh node. The bulk DMI considers the 6 NNs from every
               node an interfacial DMI is D so it only considers 4 NNs from the
               xy plane.

    OPTIONAL ARGUMENTS: -------------------------------------------------------

    dmi_type        :: 'bulk' or 'interfacial'
    name            :: Interaction name

    """

    def __init__(self, D, name='DMI', dmi_type='bulk'):
        """
        """
        self.D = D
        self.name = name
        self.jac = True
        self.dmi_type = dmi_type

    def setup(self, mesh, spin, Ms):
        super(DMI, self).setup(mesh, spin, Ms)

        # We will allow to completely specify the DMI vectors according to the
        # NNs of every lattice site, thus we need a matrix of 6 * n entries
        self.Ds = np.zeros(6 * self.n, dtype=np.float)

        if isinstance(self.D, np.ndarray) and len(self.D) == 6 * self.n:
            self.Ds = self.D
        # If we do not pass a (6 * n) array, we just create a scalar field as
        # usual and then repeat the entries 6 times so every neighbour will
        # have the same DMI per lattice site (can vary spatially)
        else:
            D_array = helper.init_scalar(self.D, self.mesh)
            self.Ds = np.repeat(D_array, 6)

        # This is from the original code:
        # self.Ds[:] = helper.init_scalar(self.D, self.mesh)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin

        if self.dmi_type == 'bulk':
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

        elif self.dmi_type == 'interfacial':
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
                "Unsppourted dmi type:{}, avaiable type: 'bulk','interfacial'.".format(self.dmi_type))

        return self.field
