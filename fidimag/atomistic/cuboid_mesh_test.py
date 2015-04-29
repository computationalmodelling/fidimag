import numpy as np
from cuboid_mesh import CuboidMesh


def allclose(a, b):
    """ close to machine precision """
    return np.allclose(a, b, rtol=1e-14, atol=1e-14)


def test_coordinates_x():
    """
       .+---+---+
     .'  .'   .'|
    +---+---+'  |        Dimensions 1 x 1 x 1.
    |   |   |   |        Cells 2 x 1 x 1.
    |   |   |   +
    |   |   | .'
    +---+---+'

    """
    mesh = CuboidMesh(1, 1, 1, 2, 1, 1)
    assert allclose(mesh.coordinates,
                    np.array(((0.25, 0.5, 0.5), (0.75, 0.5, 0.5))))


def test_coordinates_y():
    """
       .+-------+
     .+-------+'|
    +-------+'| |        Dimensions 1 x 1 x 1.
    |       | | |        Cells 1 x 2 x 1.
    |       | | +
    |       | +'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 2, 1)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.25, 0.5), (0.5, 0.75, 0.5))))


def test_coordinates_z():
    """
       .+-------+
     .'       .'|
    +-------+'  +        Dimensions 1 x 1 x 1.
    |       | .'|        Cells 1 x 1 x 2.
    +-------+'  +
    |       | .'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 1, 2)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.5, 0.25), (0.5, 0.5, 0.75))))


def test_neighbours_x():
    """
       .+---+---+
     .'  .'   .'|
    +---+---+'  |
    |   |   |   |
    | 0 | 1 |   +
    |   |   | .'
    +---+---+'

    """
    mesh = CuboidMesh(1, 1, 1, 2, 1, 1)
    assert mesh.neighbours == [{1}, {0}]


def test_neighbours_y():
    """
       .+-------+
     .+-------+'|
    +-------+'| |
    |       | |1|  
    |       |0| +
    |       | +'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 2, 1)
    assert mesh.neighbours == [{1}, {0}]


def test_neighbours_z():
    """
       .+---0---+
     .'       .'|
    +-------+'  +
    |   1   | .'|
    +-------+'  +
    |   0   | .'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 1, 2)
    assert mesh.neighbours == [{1}, {0}]


def test_neighbours_multiple():
    """
           .+-------+-------+
         .'  6    .'  7   .'|
       .+-------+-------+'  |
     .'       .'      .'|   |
    +-------+-------+'  |   +
    |       |       |   | .'|
    |   4   |   5   |   +'  |
    |       |       | .'| 3 |
    +-------+-------+'  |   +
    |       |       |   | .'
    |   0   |   1   |   +'
    |       |       | .'
    +-------+-------+'

    """
    mesh = CuboidMesh(2, 2, 2, 2, 2, 2)
    print mesh.neighbours
    assert mesh.neighbours == [{1, 2, 4}, {0, 3, 5},  # for cells 0, 1
                               {0, 3, 6}, {1, 2, 7},  # for cells 2, 3
                               {0, 5, 6}, {1, 4, 7},  # for cells 4, 5
                               {4, 2, 7}, {5, 6, 3}]  # for cells 6, 7
