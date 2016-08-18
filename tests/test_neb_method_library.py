from __future__ import print_function

"""

Tests for the NEB Method implementation

"""

import fidimag.common.neb_method as nebm
import numpy as np


def test_cartesian2spherical():
    y_cartesian = np.array([0, 0, 1,
                            0, 1, 0,
                            0, 1, 1,
                            1, 0, 0
                            ])
    y_spherical = nebm.cartesian2spherical(y_cartesian)
    theta, phi = y_spherical[::2], y_spherical[1::2]

    assert np.abs(theta[0]) < 1e-8
    assert np.abs(phi[0]) < 1e-8
    assert np.abs(theta[1] - np.pi * 0.5) < 1e-8
    assert np.abs(phi[1] - np.pi * 0.5) < 1e-8
    assert np.abs(theta[2] - np.pi * 0.25) < 1e-8
    assert np.abs(phi[2] - np.pi * 0.5) < 1e-8
    assert np.abs(theta[3] - np.pi * 0.5) < 1e-8
    assert np.abs(phi[3]) < 1e-8

    print(theta)
    print(phi)


def test_spherical2cartesian():
    A_spherical = np.array([0, 0,                 # (0, 0, 1)
                            np.pi * 0.5, np.pi,   # (-1, 0, 0)
                            np.pi * 0.25, 0,      # (1/sqrt(2), 0, 1/sqrt(2))
                            np.pi, np.pi * 0.5,   # (0, 0, -1)
                            np.pi * 0.5, 3 * np.pi * 0.5,   # (0, -1, 0)
                            ])

    A_cartesian = nebm.spherical2cartesian(A_spherical)
    A_x, A_y, A_z = A_cartesian[::3], A_cartesian[1::3], A_cartesian[2::3]

    print('A_x', A_x)
    print('A_y', A_y)
    print('A_z', A_z)

    assert np.abs(A_x[0]) < 1e-8
    assert np.abs(A_y[0]) < 1e-8
    assert np.abs(A_z[0] - 1) < 1e-8

    assert np.abs(A_x[1] - (-1)) < 1e-8
    assert np.abs(A_y[1]) < 1e-8
    assert np.abs(A_z[1]) < 1e-8

    assert np.abs(A_x[2] - 1 / np.sqrt(2)) < 1e-8
    assert np.abs(A_y[2]) < 1e-8
    assert np.abs(A_z[2] - 1 / np.sqrt(2)) < 1e-8

    assert np.abs(A_x[3]) < 1e-8
    assert np.abs(A_y[3]) < 1e-8
    assert np.abs(A_z[3] - (-1)) < 1e-8

    assert np.abs(A_x[4]) < 1e-8
    assert np.abs(A_y[4] - (-1)) < 1e-8
    assert np.abs(A_z[4]) < 1e-8


def test_linearinterpolation():
    """
    Test the linea rinterpolation of 3 pairs of spherical coordinates

                        RANGE                       INTERPOLATION
    coordinate_0  --> theta from PI/2 to PI      [PI/2 + PI/8, PI/2 + 2PI/8, PI/2 + 3PI/8]
                      phi   from 0 to 0          [ 0 0 0 ]

    coordinate_1  --> theta from PI/2 to PI/2    [0 0 0]
                      phi from PI/2 to 3PI/2     [PI/2 + PI/4, PI, PI/2 + 3PI/4]

    coordinate_2  --> We pin this one, so it should stay equal to
                      y_initial = PI, 0

    """

    y_initial = np.array([np.pi * 0.5, 0,              # coordinate 0
                          np.pi * 0.5, np.pi * 0.5,    # coordinate 1
                          np.pi * 0.5, 0               # coordinate 2
                          ])

    y_final = np.array([np.pi, 0,
                        np.pi * 0.5, 3 * np.pi * 0.5,
                        np.pi * 0.5, np.pi
                        ])

    # Pin coordinate-2
    pins = np.array([0, 0, 1])

    # Interpolate the coordinates
    interps = nebm.linear_interpolation(y_initial, y_final, 3, pins=pins)

    # Expected angles:
    expected_theta_0 = np.array([np.pi * 0.5 + np.pi * 0.125,
                                 np.pi * 0.5 + np.pi * 0.25,
                                 np.pi * 0.5 + np.pi * 0.375
                                 ])

    expected_phi_1 = np.array([np.pi * 0.5 + np.pi * 0.25,
                               np.pi * 0.5 + np.pi * 0.5,
                               np.pi * 0.5 + np.pi * 0.75
                               ])

    # Test interpolation:
    for i, theta in enumerate(interps[:, 0]):
        assert np.abs(theta - expected_theta_0[i]) < 1e-8

    for i, phi in enumerate(interps[:, 3]):
        assert np.abs(phi - expected_phi_1[i]) < 1e-8

    # The coordinates-2 is pinend so it should not change
    for i, phi in enumerate(interps[:, 5]):
        assert np.abs(phi - y_initial[5]) < 1e-8

if __name__ == "__main__":

    test_cartesian2spherical()
    test_spherical2cartesian()
    test_linearinterpolation()
