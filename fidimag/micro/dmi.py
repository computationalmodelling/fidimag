import fidimag.extensions.micro_clib as micro_clib
import numpy as np
from .energy import Energy
from fidimag.common.constant import mu_0
import fidimag.common.helper as helper
# import gc


class DMI(Energy):

    """

    Compute the Dzyaloshinskii-Moriya interaction in the micromagnetic
    framework. Currently, there are supported the following types of DMI:

        bulk    :: The energy density associated to this DMI type is:

                        w = D * M \cdot ( \nabla \times M)

                   which is found in B20 compunds or materials with
                   crystallographic class T.  Using a finite differences
                   discretisation, this term turns into an expression similar
                   to the atomistic DMI with DMI vector D_ij = -D r_ij, where
                   r_ij is the vector from the i-th mesh site to the j-th
                   neighbour

        interfacial :: The energy density of this DMI is, according to the
                       convention of Rohart et al [PRB 88, 184422 (2013)]

                            w = D * ( L_{xz}^{(x)} + L_{yz}^{(y)} )

                       where L are Lifshitz invariants. This DMI is found in
                       systems with their interface in contact with a heavy
                       metal (larger spin orbit coupling). A finite differences
                       discretisation turns this term into the equivalent
                       atomistic interfacial DMI term (with a different sign).
                       Since this DMI type is defined for interfacial systems,
                       the interaction is only defined with respect to
                       neighbouring sites in the xy plane and not between
                       layers in the z direction.

    ARGUMENTS: ----------------------------------------------------------------

    D       :: DMI vector norm which can be specified as an int, float, (X * n)
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

        if self.dmi_type == 'bulk':
            self.NN = 6
        elif self.dmi_type == 'interfacial':
            self.NN = 4

    def setup(self, mesh, spin, Ms):
        super(DMI, self).setup(mesh, spin, Ms)

        # We will allow to completely specify the DMI vectors according to the
        # NNs of every lattice site, thus we need a matrix of NN * n entries
        self.Ds = np.zeros(self.NN * self.n, dtype=np.float)

        if isinstance(self.D, np.ndarray) and len(self.D) == self.NN * self.n:
            self.Ds = self.D.astype('float')
        # If we do not pass a (NN * n) array, we just create a scalar field as
        # usual and then repeat the entries NN times so every neighbour will
        # have the same DMI per lattice site (can vary spatially)
        else:
            D_array = helper.init_scalar(self.D, self.mesh)
            self.Ds = np.repeat(D_array, self.NN)

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
                "Unsupported DMI type: {}, avaiable options: 'bulk', 'interfacial'.".format(self.dmi_type))

        return self.field
