import numpy as np


def cartesian2spherical(y_cartesian):
    """
    Convert an image of an energy band from Cartesian into spherical
    coordinates.

    Parameters
    ----------
    y_cartesian
        An image of a band of ``P + 1`` spins in Cartesian coordinates, i.e.
        with ``3 * (P + 1)`` degrees of freedom::

            [y_x0 y_y0 y_z0 y_x1 y_y1 ... y_zP]

    Returns
    -------
    numpy.ndarray
        The same image in spherical coordinates, i.e. with ``2 * (P + 1)``
        degrees of freedom::

            [y_theta0 y_phi0 y_theta1 y_phi1 ... y_phiP]

        where ``theta`` is the polar angle, ranging from 0 to PI, and
        ``phi`` is the azimuthal angle, ranging from 0 to 2 PI.
    """
    theta_phi = np.zeros((len(y_cartesian.reshape(-1, 3)), 2))

    # r = sqrt (m_x ** 2 + m_y ** 2)
    r = np.sqrt(y_cartesian[::3] ** 2 + y_cartesian[1::3] ** 2)

    # Only works if rho = sqrt(x**2 + y**2 + z**2) = 1
    # theta_phi[:, 0] = np.arccos(y_cartesian[2::3])  # theta

    theta_phi[:, 0] = np.arctan2(r, y_cartesian[2::3])  # theta
    theta_phi[:, 1] = np.arctan2(y_cartesian[1::3],
                                 y_cartesian[::3]
                                 )                      # phi

    return theta_phi.reshape(-1)


def spherical2cartesian(y_spherical):
    """
    Convert an image of an energy band from spherical into Cartesian
    coordinates.

    Parameters
    ----------
    y_spherical
        An image of a band of ``P + 1`` spins in spherical coordinates, i.e.
        with ``2 * (P + 1)`` degrees of freedom::

            [y_theta0 y_phi0 y_theta1 y_phi1 ... y_phiP]

        where ``theta`` is the polar angle, ranging from 0 to PI, and
        ``phi`` is the azimuthal angle, ranging from 0 to 2 PI.

    Returns
    -------
    numpy.ndarray
        The same image in Cartesian coordinates, i.e. with ``3 * (P + 1)``
        degrees of freedom::

            [y_x0 y_y0 y_z0 y_x1 y_y1 ... y_zP]
    """
    y_cartesian = np.zeros((len(y_spherical.reshape(-1, 2)), 3))

    theta, phi = y_spherical[::2], y_spherical[1::2]
    y_cartesian[:, 0] = np.sin(theta) * np.cos(phi)
    y_cartesian[:, 1] = np.sin(theta) * np.sin(phi)
    y_cartesian[:, 2] = np.cos(theta)

    return y_cartesian.reshape(-1)


def compute_norm(A, scale=None):
    """
    Compute the norm of the *A* array.

    Parameters
    ----------
    A
        An array with spin directions in spherical or Cartesian
        coordinates, e.g.::

            [A_theta0 A_phi0 A_theta1 A_phi1 ... A_thetaN A_phiN]

    scale
        If set to a truthy value, divide the norm by the length of *A*.

    Returns
    -------
    float
        The norm of *A*, scaled by the array size if requested.
    """

    y = np.linalg.norm(A)

    if scale:
        y = y / len(A)

    return y


def linear_interpolation_spherical(y_initial, y_final, n, pins=None):
    """
    Interpolate linearly, in spherical coordinates, between two images of an
    energy band.

    Parameters
    ----------
    y_initial
        The first image, in spherical coordinates, with the structure::

            [theta0 phi0 theta1 phi1 ...]

    y_final
        The last image, in spherical coordinates, with the same structure as
        *y_initial*.
    n
        Number of interpolations to generate between the two images.
    pins
        An array or list of 0s and 1s, representing unpinned and pinned
        sites respectively, whose entries are not interpolated. Since there
        is one entry per site, the *pins* array must have HALF the length of
        the *y* arrays::

            [pin0 pin1 ...]

        If ``None``, no site is pinned.

    Returns
    -------
    numpy.ndarray
        An array of shape ``(n, len(y_initial))``, where every row is one of
        the interpolated images.
    """

    # We will generate n copies of the y_initial array, using rows
    # For this, we use Numpy's broadcasting. For example,
    # if y_initial=[1 3 4 5], then:
    #
    #        [ [0]     + [1 3 4 5]  = [ [1  3  4  5]
    #          [0] ]                    [1  3  4  5] ]
    interpolations = np.zeros((n, 1))
    interpolations = interpolations + y_initial

    # We will not interpolate pinned spins
    if pins is None:
        #  Just use half the length of y_initial
        pins = np.zeros(len(y_initial[::2]))

    # Since we have a pin index per every PAIR of coordinates, we copy every
    # entry. For example: [1 0] --> [1 1 0 0]
    # and we change only unpinned spins (0)
    _filter = np.repeat(pins, 2) == 0

    # y_initial_spherical = cartesian2spherical(y_initial)
    # y_final_spherical = cartesian2spherical(y_final)

    # dy_spherical = ((y_final_spherical - y_initial_spherical) / (n + 1))
    dy = (y_final - y_initial) / (n + 1)

    for i in range(1, n + 1):
        interpolations[i - 1][_filter] = (y_initial + i * dy)[_filter]

    return interpolations


def interpolation_Rodrigues_rotation(y_initial, y_final, n, pins=None):
    """
    Interpolate between two images of an energy band using the Rodrigues
    rotation formula, in order to generate a sequence of NEBM images.

    The angles of the interpolations are the angles between corresponding
    spins of the *y_initial* and *y_final* vectors, computed as the arccos
    of their dot products.

    Parameters
    ----------
    y_initial
        The first image, in Cartesian coordinates, with the structure::

            [mx1 my1 mz1, mx2 my2 ...]

    y_final
        The last image, in Cartesian coordinates, with the same structure as
        *y_initial*.
    n
        Number of interpolations to generate between the two images.
    pins
        An array or list of 0s and 1s, representing unpinned and pinned
        sites respectively, whose entries are not interpolated. Since there
        is one entry per site, the *pins* array must have a THIRD of the
        length of the *y* arrays::

            [pin0 pin1 ...]

        If ``None``, no site is pinned.

    Returns
    -------
    numpy.ndarray
        An array of shape ``(n, len(y_initial))``, where every row is one of
        the interpolated images.

    References
    ----------
    Bessarab et al., Computer Physics Communications 196 (2015) 335-347
    """

    # We will generate n copies of the y_initial array, using rows
    # For this, we use Numpy's broadcasting. For example,
    # if y_initial=[1 3 4 5], then:
    #
    #        [ [0]     + [1 3 4 5]  = [ [1  3  4  5]
    #          [0] ]                    [1  3  4  5] ]
    interpolations = np.zeros((n, 1))
    interpolations = interpolations + y_initial

    # We will not interpolate pinned spins ------------------------------------
    if pins is None:
        # Only use 1/3 of the length of y_initial (1 pin per mesh/lattice site)
        pins = np.zeros(len(y_initial[::3]))

    # Since we have a pin index per every TRIAD of coordinates, we copy every
    # entry. For example: [1 0] --> [1 1 1 0 0 0]
    # and we change only unpinned spins (0)
    _filter = np.repeat(pins, 3) == 0

    # -------------------------------------------------------------------------

    # We will perform the calculations for every spin at once
    y_initial.shape = (-1, 3)
    y_final.shape = (-1, 3)

    # The cross products of corresponding spins in the initial and final images
    yi_cross_yf = np.cross(y_initial, y_final)
    # This should only be an array of ones:
    yi_cross_yf_norm = np.linalg.norm(yi_cross_yf, axis=1)

    # The rotation axis is just the normalised cross product defined before
    rot_axis = np.zeros_like(yi_cross_yf)
    fltr = yi_cross_yf_norm > 0
    rot_axis[fltr] = yi_cross_yf[fltr] / yi_cross_yf_norm[fltr][:, np.newaxis]
    rot_axis = np.cross(rot_axis, y_initial)

    # The angles between corresponding spins
    yi_yf_angle = np.arccos(np.sum(y_initial * y_final, axis=1))

    for i in range(1, n + 1):
        dangle = i * yi_yf_angle / (n + 1)

        # Rodrigues formulae for the i-th interpolation
        interp = (y_initial * np.cos(dangle)[:, np.newaxis]
                  + rot_axis * np.sin(dangle)[:, np.newaxis]
                  ).reshape(-1)
        interpolations[i - 1][_filter] = interp[_filter]

    y_initial.shape = (-1)
    y_final.shape = (-1)

    return interpolations


def m_to_zero_nomaterial(image_cartesian, sim):
    """
    Set the spin direction of the sites with no material to ``[0, 0, 0]``.

    Parameters
    ----------
    image_cartesian
        An image of a band in Cartesian coordinates::

            [mx_0 my_0 mz_0 mx_1 my_1 ... mz_(P-1)]

    sim
        A fidimag simulation object, from which the ``Ms`` or ``mu_s``
        values are taken to filter the spin directions.

    Returns
    -------
    numpy.ndarray
        A copy of the image with the sites with no material set to zero.
    """
    image_reshape = np.copy(image_cartesian.reshape(-1, 3))
    _filter = sim._magnetisation == 0
    image_reshape[_filter] = np.array([0., 0., 0.])
    return image_reshape.reshape(-1)
