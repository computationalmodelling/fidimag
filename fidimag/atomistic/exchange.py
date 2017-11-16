import fidimag.extensions.clib as clib
import numpy as np
from .energy import Energy


class Exchange(Energy):

    """
    This class provides the Exchange Interaction energy term, defined as

                  __        ->      ->
         E =  -  \    J_ij  S_i  *  S_j
                 /__
                <i, j>
                i != j

    where J_ij is the exchange tensor, S_i and S_j are the total spin vectors
    at the i-th and j-th lattice sites, and <i, j> means counting the
    interaction between neighbouring spins only once (notice that there is no
    factor of 2 associated with J)

    In general, J_ij is space dependent and must be specified for every
    neighbour.  For a homogeneous material, J_ij is a diagonal tensor with
    constant magnitude, i.e. J_ij -> J.


    OPTIONAL ARGUMENTS: -------------------------------------------------------

        J               :: The exchange tensor which can be  a number (same
                           magnitude for every neighbour at every lattice site)
                           or a space dependent function that returns 6
                           components, one for every nearest neighbour (NN).
                           For a square lattice the NNs are defined in 3D as:
                           [-x +x -y +y -z +z], thus the exchange components
                           are specified in that order. In a hexagonal lattice
                           the NNs are only in a 2D plane as the cardinal
                           positions: [W E NE SW NW SE].

        name            :: Interaction name

    USAGE: --------------------------------------------------------------------

    For a homogeneous material, it can be specified in a simulation object
    *Sim* as

            Sim.add(Exchange(J))

    where J is a float.

    Otherwise, the exchange tensor is spatial dependent, and must be specified
    using a function. For a cubic mesh, for example, if we want a linear
    dependence only with in plane components:

            def my_exchange(pos):
                J = 1 * meV
                x, y, z = pos[0], pos[1], pos[2]

                return (x * J, x * J, y * J, y * J, 0, 0)

            Sim.add(Exchange(my_exchange))

    DEV NOTES: ----------------------------------------------------------------

    * If a float or int is passed as the Exchange constant J, this class will
    use the *compute_exchange_field* C function (see lib/exch.c), which assumes
    a uniform exchange. This C function does not take J as an array thus it
    will not call array elements to compute the neighbours contribution but
    it will only use a constant, thus it should be faster

    """

    def __init__(self, J, name='Exchange'):
        self.J = J
        self.name = name
        self.jac = False

    def setup(self, mesh, spin, mu_s):
        super(Exchange, self).setup(mesh, spin, mu_s)

        if isinstance(self.J, (int, float)):
            self.Jx = float(self.J)
            self.Jy = float(self.J)
            self.Jz = float(self.J)

            self.compute_field = self.compute_field_uniform

        # TODO: Add option to pass numpy arrays
        elif hasattr(self.J, '__call__'):
            self._J = np.zeros(self.neighbours.shape)
            n = self.mesh.n
            for i in range(n):
                value = self.J(self.coordinates[i])
                if len(value) == 6:
                    self._J[i, :] = value[:]
                else:
                    raise Exception('The given spatial function for J is not acceptable!')
            pass

            self.compute_field = self.compute_field_spatial

        elif isinstance(self.J, list) and len(self.J) == self.mesh.shells:
            self._J = [0. for i in range(8)]
            for i in range(len(self.J)):
                self._J[i] = float(self.J[i])

            self.compute_field = self.compute_field_full()

    def compute_field_spatial(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_exchange_field_spatial(m,
                                            self.field,
                                            self.energy,
                                            self._J,
                                            self.neighbours,
                                            self.n,
                                            self.n_ngbs
                                            )

        return self.field * self.mu_s_inv

    def compute_field_uniform(self, t=0, spin=None):

        if spin is not None:
            m = spin
        else:
            m = self.spin

        clib.compute_exchange_field(m,
                                    self.field,
                                    self.energy,
                                    self.Jx,
                                    self.Jy,
                                    self.Jz,
                                    self.neighbours,
                                    self.n,
                                    self.n_ngbs
                                    )

        return self.field * self.mu_s_inv

    def compute_field_full(self):
        pass


class UniformExchange(Exchange):

    """
    For compatibility we leave this class which was merged into the
    Exchange class
    """

    def __init__(self, J=0, name='UniformExchange'):
        super(UniformExchange, self).__init__(J, name=name)
