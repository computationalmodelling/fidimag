import fidimag.extensions.clib
import fidimag.extensions.micro_clib
import numpy as np


def skyrmion_number_at_centre(sim):
    """
    Returns the skyrmion number calculated on the x-y plane at the central z
    co-ordinate.

    The current fidimag skyrmion number function finds the skyrmion number for
    the first z co-ordinate in index order (unintentionally?). So the approach
    this function uses is to slice the spin array in the simulation object to
    only include the x-y plane as previously described.

    Arguments:
      sim: LLG object with cuboidal mesh.

    Returns:
      skyrmionNumber (a float)
    """

    # Find the "length" of a slice in terms of the spin array.
    xyLength = sim.mesh.nx * sim.mesh.ny * 3  # The 3 represents (x,y,z).

    # Find the layer number that corresponds to the centre of the mesh,
    # preferring the lower layer in a tie.
    zCentre = sim.mesh.nz / 2.
    if sim.mesh.nz % 2 == 0:
        zCentre -= 0.5
    zCentre -= 0.5

    # Obtain slice of spins cleverly. This also works if the domain is flat.
    spinSlice = sim.spin[int(xyLength * zCentre):int(xyLength * (zCentre + 1))]

    # Compute the skyrmion number for our spin slice instead.
    if sim._micromagnetic is True:
        return fidimag.extensions.micro_clib.compute_skyrmion_number(\
               spinSlice, sim._skx_number, sim.mesh.nx, sim.mesh.ny,
               sim.mesh.nz, sim.mesh.neighbours)
    else:
        return fidimag.extensions.clib.compute_skyrmion_number(\
               spinSlice, sim._skx_number, sim.mesh.nx, sim.mesh.ny,
               sim.mesh.nz, sim.mesh.neighbours)


def skyrmion_number_lee(sim):
    """
    Returns the skyrmion number calculated from a 3D sample, as defined in:

    Lee, M, Kang, W, Onose, Y, et. al (2009) "Unusual Hall effect anomaly in
    MnSi under pressure". PRL.

    Arguments:
      sim: LLG object with cuboidal mesh.

    Returns:
      skyrmionNumber (a float)
    """

    # Find the "length" of a slice in terms of the spin array.
    xyLength = sim.mesh.nx * sim.mesh.ny * 3  # The 3 represents (x, y, z).

    # Create an array to store skyrmion number values for each z slice.
    skyrmionNumbers = np.ndarray(sim.mesh.nz)

    # For each slice, compute the skyrmion number.
    for zI in range(len(skyrmionNumbers)):
        spinSlice = sim.spin[xyLength * zI:xyLength * (zI + 1)]

        if sim._micromagnetic is True:
            skyrmionNumbers[zI] = fidimag.extensions.micro_clib.compute_skyrmion_number(\
                                  spinSlice, sim._skx_number, sim.mesh.nx,
                                  sim.mesh.ny, sim.mesh.nz,
                                  sim.mesh.neighbours)
        else:
            skyrmionNumbers[zI] = fidimag.extensions.clib.compute_skyrmion_number(\
                                  spinSlice, sim._skx_number, sim.mesh.nx,
                                  sim.mesh.ny, sim.mesh.nz,
                                  sim.mesh.neighbours)

    # Return the average. This is equivalent to the integral in the equation.
    return skyrmionNumbers.mean()
