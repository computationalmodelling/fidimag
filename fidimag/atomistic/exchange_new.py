from field import scalar_field, vector_field

class Exchange(object):
    """
    Uniform exchange interaction.

    """
    def __init__(self, mesh, J, name="exchange"):
        self.mesh = mesh
        # J is a scalar field, because that's good enough for uniform
        # exchange. Could be changed trivially by making J a vector field
        # and changing the J = self.J[c_i] line in the compute method.
        self.J = scalar_field(mesh, J)
        # field and energy are just np.arrays. In principle, they could be
        # created in the scope of the compute method, but we want to avoid
        # re-creating the objects in memory unnecessarily.
        self.field = vector_field(mesh, 0)
        self.energy = scalar_field(mesh, 0)
        self.in_jacobian = True

    def compute(self, spins):
        """
        Compute the exchange field and the exchange energy.

        """
        for c_i in mesh.cells():
            #TODO: cythonise this loop fiercely
            # like shown at http://docs.cython.org/src/tutorial/numpy.html
            fx = 0
            fy = 0
            fz = 0
            J = self.J[c_i]
            for c_j in mesh.neighbours[c_i]:
                fx += J * spins[c_j, 0]
                fy += J * spins[c_j, 1]
                fz += J * spins[c_j, 2]
            self.field[c_i] = (fx, fy, fz)
            self.energy[c_i] = -0.5 * (  fx * spins[c_i, 0]
                                       + fy * spins[c_i, 1]
                                       + fz * spins[c_i, 2])
