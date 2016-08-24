from __future__ import print_function
from __future__ import division
import numpy as np


def cartesian2spherical(y_cartesian):
    """
    y_cartesian     :: [y_x0 y_y0 y_z0 y_x1 y_y1 ...]
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
    y_cartesian = np.zeros((len(y_spherical.reshape(-1, 2)), 3))

    theta, phi = y_spherical[::2], y_spherical[1::2]
    y_cartesian[:, 0] = np.sin(theta) * np.cos(phi)
    y_cartesian[:, 1] = np.sin(theta) * np.sin(phi)
    y_cartesian[:, 2] = np.cos(theta)

    return y_cartesian.reshape(-1)


def compute_norm(A, scale=None):
    """

    Compute the norm of the *A* array, which contains spin directions in
    Spherical coordinates,

    A = [ A_theta0 A_phi0 A_theta1 A_phi1 ... A_thetaN A_phiN]

    If the absolute value of a component is larger than PI, we redefine the
    difference to be smaller than PI

    We scale the norm by the array size

    """

    y = np.copy(A)

    if scale:
        y = np.sqrt(np.sum(y ** 2.)) / len(y)
    else:
        y = np.sqrt(np.sum(y ** 2.))

    return y


def linear_interpolation_spherical(y_initial, y_final, n, pins=None):
    """

    This function returns a (n, len(y_initial)) array, where every row is
    an interpolation of the coordinates in y_initial to y_final. The
    interpolation is made in spherical coordinates

    ARGUMENTS:

    y_initial, y_final      :: In Spherical coordinates with the structure:

                                    [ theta0 phi0 theta1 phi1 ...]

    OPTIONAL:

    pins                    :: An array or list with 0s and 1s,
                               representing unpinned/pinned coordinates of
                               any of the *y* arrays, respectively. Thus,
                               the *pins* array must have HALF the length
                               of *y*:

                                    [pin0 pin1  ...  ]

    """

    # We will generate n copies of the y_initial array, using rows
    # For this, we use Numpy's broadcasting. For example,
    # if y_initial=[1, 3 4, 5], then:
    #
    #        [ [0]     + [1 3 4 5]  = [ [1  3  4  5]
    #          [0] ]                    [1  3  4  5] ]
    interpolations = np.zeros((n, 1))
    interpolations = interpolations + y_initial

    # We will not interpolate pinned spins
    if pins is None:
        #  Just use half the length of y_initial
        pins = np.zeros(len(y_initial[::2]))

    # Since we have a pin index per every PAIR of coordinates, we copy very
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
