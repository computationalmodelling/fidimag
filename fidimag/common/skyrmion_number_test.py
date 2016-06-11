from fidimag.common.skyrmion_number import skyrmion_number_at_centre
from fidimag.common.skyrmion_number import skyrmion_number_lee
import copy
import fidimag
import numpy as np


def skyrmion_centre_z(mesh):
    """
    Produces a function defining a magnetisation texture with a skyrmion on the
    centre layer of the mesh only, and ferromagnetic otherwise.

    Arguments:
        mesh: Cuboid mesh object from Fidimag.

    Returns a function object.
    """

    def m_skyrmion_unrealistic(pos):
        """
        Produces a skyrmion on the centre layer of the mesh only, and
        ferromagnetic otherwise.

        Arguments:
            pos: 3-element container of floats denoting the position to probe
                 the magnetisation of the skyrmion.

        Returns a three-element numpy array representing the magnetisation of
        the skyrmion.
        """

        # Find the layer number that corresponds to the centre of the mesh,
        # preferring the lower layer in a tie.
        zCentre = mesh.nz / 2.
        if mesh.nz % 2 == 0:
            zCentre -= 0.5

        # Find the co-ordinate of that layer.
        zCoord = mesh.z0 + mesh.dz * zCentre

        # Logic as described in docstring.
        if abs(pos[2] - zCoord) < 1e-6:

            # This skyrmion is centred in the centre of the mesh, and has
            # radius big enough to fill the mesh, assuming the x and y lengths
            # are the same.
            skyrmionCentre = [mesh.x0 + mesh.dx * mesh.nx / 2.,
                              mesh.y0 + mesh.dy * mesh.ny / 2.]
            skyrmionRadius = min(mesh.dx * mesh.nx / 2.,
                                 mesh.dy * mesh.ny / 2.)
            return m_skyrmion_centre(pos, skyrmionCentre, skyrmionRadius)

        else:
            return np.array([0., 0., -1.], dtype="float64")

    # Return the function object, as promised.
    return m_skyrmion_unrealistic


def skyrmion_trouser_leg(mesh):
    """
    Produces a magnetisation texture with a skyrmion invariant in thickness.

    Arguments:
        mesh: Cuboid mesh object from Fidimag.

    Returns a function object.
    """

    # This skyrmion is centred in the centre of the mesh, and has radius big
    # enough to fill the mesh, assuming the x and y lengths are the same.
    skyrmionCentre = [mesh.x0 + mesh.dx * mesh.nx / 2.,
                      mesh.y0 + mesh.dy * mesh.ny / 2.]
    skyrmionRadius = min(mesh.dx * mesh.nx / 2., mesh.dy * mesh.ny / 2.)

    def m_skyrmion_trouser_leg(pos):
        """
        Produces a skyrmion invariant in thickness.

        Arguments:
            pos: 3-element container of floats denoting the position to probe
                 the magnetisation of the skyrmion.

        Returns a three-element numpy array representing the magnetisation of
        the skyrmion.
        """
        return m_skyrmion_centre(pos, skyrmionCentre, skyrmionRadius)

    # Return the function object, as promised.
    return m_skyrmion_trouser_leg


def m_skyrmion_centre(pos, skyrmionCentre, skyrmionRadius):
    """
    Returns the magnetisation of a skyrmion at a given position.

    Arguments:
        pos: 3-element container of floats denoting the position to probe the
             magnetisation of the skyrmion.
        skyrmionCentre: 2-element list of floats denoting the centre-point
                        of the skyrmion.
        skyrmionRadius: Float denoting the radius of the skyrmion.

    Returns a three-element numpy array representing the magnetisation of the
    skyrmion.
    """

    # Convert co-ordinate to polar form centred around the skyrmion.
    loc = copy.deepcopy(pos[:2])
    loc -= np.array(skyrmionCentre)
    r = np.linalg.norm(loc)

    # Check to see if this vector is within this circle. If it isn't,
    # check the next centre by skipping the rest of this iteration.
    if r > skyrmionRadius:
        return np.array([0., 0., -1.], dtype="float64")

    # Convert position into cylindrical form, using r defined
    # previously, and "t" as the argument.
    if abs(loc[0]) < 1e-6:
        if abs(loc[1]) < 1e-6:
            t = 0
        elif loc[1] > 0:
            t = np.pi / 2.
        else:
            t = -np.pi / 2.
    else:
        t = np.arctan2(loc[1], loc[0])

    # Define vector components inside the skyrmion:
    mz = np.cos(np.pi * r / skyrmionRadius)
    mt = -np.sin(np.pi * r / skyrmionRadius)

    # Convert to cartesian form and normalize.
    mx = -np.sin(t) * mt
    my = np.cos(t) * mt
    out = np.array([mx, my, mz], dtype="float64")
    return out / np.linalg.norm(out)


# Let the tests begin!
class Case():
    def __init__(self):

        # Building mesh and simulation objects.
        layers = 3
        self.mesh = fidimag.common.CuboidMesh(dx=1, dy=1, dz=1,
                                              nx=20, ny=20, nz=layers,
                                              x0=0, y0=0, z0=0,
                                              unit_length=1e-9)
        self.sim = fidimag.micro.Sim(self.mesh)

        # Absolute tolerance used for comparisons.
        self.aTol = 5e-2


def test_skyrmion_number_monolayer():
    """
    Tests skyrmion number function calculations with a mesh initialised uniform
    everywhere, except a skyrmion in the centre z-layer.

    We expect that the clib skyrmion number calculation only finds the skyrmion
    number for the first z-slice of a mesh, that the centre skyrmion number
    function finds the skyrmion at the centre, and that the 3D skyrmion number
    function takes an average.

    Returns nothing.
    """

    testCase = Case()

    print("Case: skyrmion in one layer only.")
    testCase.sim.set_m(skyrmion_centre_z(testCase.mesh))

    clibSk = testCase.sim.skyrmion_number()
    print("CLib calculates the skyrmion number as: {:1.2f}.".format(clibSk))

    centreSk = skyrmion_number_at_centre(testCase.sim)
    print("I calculate the skyrmion number at the centre as: {:1.2f}."
          .format(centreSk))

    leeSk = skyrmion_number_lee(testCase.sim)
    print("I calculate the 3D skyrmion number as: {:1.2f}.".format(leeSk))

    assert abs(clibSk - 0) < testCase.aTol
    assert abs(centreSk - 1.) < testCase.aTol
    assert abs(leeSk - 1 / float(testCase.sim.mesh.nz)) < testCase.aTol


def test_skyrmion_number_multilayer():
    """
    Tests skyrmion number function calculations with a mesh initialised with a
    skyrmion throughout all z-layers of the mesh.

    We expect that the clib skyrmion number calculation only finds the skyrmion
    number for the first z-slice of a mesh, that the centre skyrmion number
    function finds the skyrmion at the centre, and that the 3D skyrmion number
    function takes an average.

    Returns nothing.
    """

    testCase = Case()

    print("Case: skyrmion consistent across layers (skyrmion trouser leg)")
    testCase.sim.set_m(skyrmion_trouser_leg(testCase.mesh))

    clibSk = testCase.sim.skyrmion_number()
    print("CLib calculates the skyrmion number as: {:1.2f}.".format(clibSk))

    centreSk = skyrmion_number_at_centre(testCase.sim)
    print("I calculate the skyrmion number at the centre as: {:1.2f}."
          .format(centreSk))

    leeSk = skyrmion_number_lee(testCase.sim)
    print("I calculate the 3D skyrmion number as: {:1.2f}.".format(leeSk))

    assert abs(clibSk - 1.) < testCase.aTol
    assert abs(centreSk - 1.) < testCase.aTol
    assert abs(leeSk - 1.) < testCase.aTol


# Run tests if this script was executed.
if __name__ == "__main__":
    test_skyrmion_number_monolayer()
    test_skyrmion_number_multilayer()
