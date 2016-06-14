import fidimag.extensions.clib
import fidimag.extensions.micro_clib
import numpy as np


def skyrmion_number_slice(sim, at=None, zIndex=None):
    """
    Returns the skyrmion number calculated on the x-y plane at a specific z
    co-ordinate.

    The current fidimag skyrmion number function finds the skyrmion number for
    the first z co-ordinate in index order (unintentionally?). So the approach
    this function uses is to slice the spin array in the simulation object to
    only include the x-y plane as previously described.

    Arguments:
      sim: LLG object with cuboidal mesh.

    Optional Arguments:
      If neither of these arguments are specified, this function raises a
      ValueError.

      at: None, "top", "centre", or "bottom" denoting where to slice in
          British English. https://en.wikipedia.org/wiki/Strike_It_Lucky
      zIndex: Integer (not a float), or None denoting the index at which to
          slice. Overrides "at" if not None. Bound by zIndex~[0, sim.mesh.nz-1]

    Returns:
      skyrmionNumber (a float)
    """

    # Deal with inputs.
    if zIndex is None:
        if at is "bottom":
            zIndex = 0

        elif at is "centre":
            # Find the layer number that corresponds to the centre of the mesh,
            # preferring the lower layer in a tie.
            zIndex = sim.mesh.nz / 2.
            if sim.mesh.nz % 2 == 0:
                zIndex -= 0.5
            zIndex -= 0.5

        elif at is "top":
            zIndex = sim.mesh.nz - 1  # Indeces start at 0.

        else:
            raise ValueError("[skyrmion_number_slice]: Either "
                             "'at=bottom|centre|top' or 'zIndex'=[integer] "
                             "must be passed as arguments.")

    elif zIndex > sim.mesh.nz - 1 or zIndex < 0:
        raise ValueError("[skyrmion_number_slice]: 'zIndex={}' must be an "
                         "integer between 0 and sim.mesh.nz - 1, which is "
                         "'{}' for this simulation object."
                         .format(zIndex, sim.mesh.nz - 1))

    # Find the "length" of a slice in terms of the spin array.
    xyLength = sim.mesh.nx * sim.mesh.ny * 3  # The 3 represents (x,y,z).

    # Obtain slice of spins cleverly. This also works if the domain is flat
    # (i.e. nz=1, zIndex=0)
    spinSlice = sim.spin[int(xyLength * zIndex):int(xyLength * (zIndex + 1))]

    # Compute the skyrmion number for our spin slice.
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
