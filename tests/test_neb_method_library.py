from __future__ import print_function

"""

Tests for the NEB Method implementation

"""

import fidimag.extensions.nebm_geodesic_clib as nebm_geodesic
import fidimag.extensions.nebm_cartesian_clib as nebm_cartesian
import fidimag.common.nebm_tools as nebm_tools
# import fidimag.common.nebm_spherical as nebm_spherical
import numpy as np


def test_cartesian2spherical():
    y_cartesian = np.array([0, 0, 1,
                            0, 1, 0,
                            0, 1, 1,
                            1, 0, 0
                            ])
    y_spherical = nebm_tools.cartesian2spherical(y_cartesian)
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

    A_cartesian = nebm_tools.spherical2cartesian(A_spherical)
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
    interps = nebm_tools.linear_interpolation_spherical(y_initial, y_final,
                                                        3, pins=pins)

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


def test_geodesic_distance():
    """
    We will measure the distance between a spin in the +z direction and a spin
    at the -z direction. The great-circle distance (distance on the unit
    sphere) should be PI, since the perimeter is 2PI

    The other distance is between a spin in the [0, 1, 0] direction and another
    one in the [-1, 0, 0] direction. The arc length should be PI/2
    """

    # Degrees of freedom in Cartesian coordinates
    A = np.array([0, 0, 1.])
    B = np.array([0, 0, -1.])
    n_dofs_image = 3

    # The material array only indicates the lattice/mesh sites where
    # there is material (mu_s or M_s larger than zero)
    material = np.array([1]).astype(np.int32)
    n_dofs_image_material = 3

    d = nebm_geodesic.geodesic_distance(A, B, n_dofs_image,
                                        material, n_dofs_image_material
                                        )
    assert np.abs(d - np.pi) < 1e-5

    # The other two spins:
    A = np.array([0, 1, 0.])
    B = np.array([-1, 0, 0.])

    d = nebm_geodesic.geodesic_distance(A, B, n_dofs_image,
                                        material, n_dofs_image_material
                                        )
    assert np.abs(d - np.pi * 0.5) < 1e-5


# Cartesian NEBM library: common/neb_method/nebm_cartesian_lib.c
def test_project_vectors_into_image():
    """
    Project 3 vectors in array *a into an image made of
    9 degrees of freedom = 3 vector
    """

    a = np.arange(9, dtype=np.float)
    image = np.array([0., 0, 1,
                      1, 0, 0,
                      0, 0, 1])

    # Project array of vectors a, which is updated in the C library
    nebm_cartesian.project_vector(a, image, 9)

    # Here we do the projection manually --------------------------------------
    # Dot products for every vector
    a_dot_image = [2., 3., 8.]
    expected_res = np.zeros(9)
    # Projections:
    expected_res[:3] = [0, 1, 2 - a_dot_image[0]]
    expected_res[3:6] = [3 - a_dot_image[1], 4, 5]
    expected_res[6:] = [6, 7, 8 - a_dot_image[2]]
    # print(expected_res)
    # print(a)
    # -------------------------------------------------------------------------

    for i in range(9):
        assert a[i] == expected_res[i]


def test_project_images_into_image():
    """

    Project the vectors from an array of images a, into images

    The C library function does NOT calculate projections of extrema images,
    i.e. image 0 and N, thus we pad two extra arrays at the end and beginning
    of the arrays

    We construct array "a" made of 4 images, each having 3 vectors
    Then we construct array "images" with 4 images
    Then we project vectors from "a" into "images" and compare results

    """

    a = np.zeros(18 + 18, dtype=np.float)
    a[9:27] = np.arange(18, dtype=np.float)
    images = np.array([0., 0, 0, 0, 0, 0, 0, 0, 0,
                       0., 0, 1, 1, 0, 0, 0, 0, 1,
                       0., 1, 0, 0, 0, 1, 1, 1, 0,
                       0., 0, 0, 0, 0, 0, 0, 0, 0,
                       ])

    # Project array of images a, which is updated in the C library
    nebm_cartesian.project_images(a, images,
                                  # num of images and num of dofs per image:
                                  4, 9)

    # Here we do the projections manually -------------------------------------
    # Dot products for every vector
    a_dot_images = [2., 3., 8.,
                    10., 14., 15. + 16.]
    expected_res = np.zeros(18 + 18, dtype=np.float)
    # Projections:
    expected_res[9:12] = [0, 1, 2 - a_dot_images[0]]
    expected_res[12:15] = [3 - a_dot_images[1], 4, 5]
    expected_res[15:18] = [6, 7, 8 - a_dot_images[2]]
    # Second image
    expected_res[18:21] = [9, 10 - a_dot_images[3], 11]
    expected_res[21:24] = [12, 13, 14 - a_dot_images[4]]
    expected_res[24:27] = [15 - a_dot_images[5], 16 - a_dot_images[5], 17]
    # print(expected_res)
    # print(a)
    # -------------------------------------------------------------------------

    for i in range(36):
        assert a[i] == expected_res[i]


if __name__ == "__main__":

    test_cartesian2spherical()
    test_spherical2cartesian()
    test_linearinterpolation()
    test_geodesic_distance()
    test_project_vectors_into_image()
    test_project_images_into_image()
