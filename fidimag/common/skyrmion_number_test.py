#!/usr/local/bin/ipython

from skyrmion_number_lib import skyrmion_number_at_centre, skyrmion_number_lee
import copy
import numpy as np
import fidimag  # Revision 9539fbc


def skyrmion_centre_z(sim):
    """
    Produces a hackish skyrmion function, with a skyrmion on the centre later
    only, and ferromagnetic otherwise.
    """

    def m_hackish(pos):
        """
        Produces a skyrmion on the centre layer only, and ferromagnetic
        otherwise.
        """

        # Find the layer number that corresponds to the centre of the mesh,
        # preferring the lower layer in a tie.
        zCentre = sim.mesh.nz / 2.
        if sim.mesh.nz % 2 == 0:
            zCentre -= 0.5

        # Find the co-ordinate of that layer.
        zCoord = sim.mesh.z0 + sim.mesh.dz * zCentre

        # Logic as described in docstring.
        if abs(pos[2] - zCoord) < 1e-6:
            return m_skyrmion_centre(pos)
        else:
            return np.array([0., 0., 1.], dtype="float64")
    return m_hackish


def skyrmion_trouser_leg(sim):
    """
    Produces a hackish skyrmion function, with a skyrmion invariant in
    thickness.
    """

    return m_skyrmion_centre


def m_skyrmion_centre(pos):
    """
    Produces vector field mimicing a skyrmion of radius 30 nanometres.
    """
    skyrmionRadius = 30.

    loc = copy.deepcopy(pos)

    # Place it in the centre of our sample.
    loc -= np.array([50, 50, 0])

    # Find the radius component in cylindrical form.
    r = np.linalg.norm(loc)

    # Check to see if this vector is within this circle. If it isn't,
    # check the next centre by skipping the rest of this iteration.
    if r > skyrmionRadius:
        return np.array([0., 0., 1.], dtype="float64")

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
    mz = -np.cos(np.pi * r / skyrmionRadius)
    mt = np.sin(np.pi * r / skyrmionRadius)

    # Convert to cartesian form and normalize.
    mx = -np.sin(t) * mt
    my = np.cos(t) * mt
    out = np.array([mx, my, mz], dtype="float64")
    return out / np.linalg.norm(out)


layers = 31  # Should be odd.
mesh = fidimag.common.CuboidMesh(dx=1, dy=1, dz=1,
                                 nx=100, ny=100, nz=layers,
                                 x0=0, y0=0, z0=0, unit_length=1e-9,
                                 periodicity=[True, True, False])

sim = fidimag.micro.Sim(mesh)
tolerance = 1e-1

# Check expected results for a skyrmion on one layer only.
print("\nFirst case: skyrmion in one layer only.")
sim.set_m(skyrmion_centre_z(sim))
sim.save_vtk()

fidimagSk = sim.skyrmion_number()
print("Fidimag calculates the skyrmion number as: {:1.2f}.".format(fidimagSk))

centreSk = skyrmion_number_at_centre(sim)
print("I calculate the skyrmion number at the centre as: {:1.2f}."
      .format(centreSk))

leeSk = skyrmion_number_lee(sim)
print("I calculate the 3d skyrmion number as: {:1.2f}.".format(leeSk))

assert abs(fidimagSk - 0) < tolerance
assert abs(centreSk - 1.) < tolerance
assert abs(leeSk - 1 / float(layers)) < tolerance


# Check expected results for a skyrmion invariant in z.
print("\nSecond case: skyrmion consistent across layers.")
sim.set_m(skyrmion_trouser_leg(sim))

fidimagSk = sim.skyrmion_number()
print("Fidimag calculates the skyrmion number as: {:1.2f}.".format(fidimagSk))
assert abs(fidimagSk - 1.) < tolerance

centreSk = skyrmion_number_at_centre(sim)
print("I calculate the skyrmion number at the centre as: {:1.2f}."
      .format(centreSk))
assert abs(centreSk - 1.) < tolerance

leeSk = skyrmion_number_lee(sim)
print("I calculate the 3d skyrmion number as: {:1.2f}.".format(leeSk))
assert abs(leeSk - 1.) < tolerance
