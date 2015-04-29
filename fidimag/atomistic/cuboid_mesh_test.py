import numpy as np
from cuboid_mesh import \
        CuboidMesh, \
        vector_from_coordinates, \
        scalar_from_coordinates


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


def test_iterate_over_cells():
    mesh = CuboidMesh(1, 1, 1, 2, 2, 2)
    for c_i in mesh.cells():
        print "This is cell #{}, I have neighbours {}.".format(c_i, mesh.neighbours[c_i])


def test_iterate_over_cells_and_neighbours():
    mesh = CuboidMesh(1, 1, 1, 2, 2, 2)
    for c_i in mesh.cells():
        print "I am cell #{}.".format(c_i)
        for c_j in mesh.neighbours[c_i]:
            print "\tAnd I am its neighbour, cell #{}!".format(c_j)


def test_scalar_field_shape():
    mesh = CuboidMesh(1, 1, 1, 2, 3, 4)
    expected_nb_cells = 2 * 3 * 4
    expected_shape_for_scalar_field = (expected_nb_cells, 1)
    assert mesh.scalar_shape() == expected_shape_for_scalar_field
    f = np.zeros(mesh.scalar_shape())  # usage example


def test_vector_field_shape():
    mesh = CuboidMesh(1, 1, 1, 2, 3, 4)
    expected_nb_cells = 2 * 3 * 4
    expected_shape_for_vector_field = (expected_nb_cells, 3)
    assert mesh.vector_shape() == expected_shape_for_vector_field
    m = np.zeros(mesh.vector_shape())  # usage example


def test_initialise_vector():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    v = vector_from_coordinates(lambda r: 2 * r, mesh)
    assert allclose(v, np.array((1, 1, 1)))


def test_initialise_scalar():
    mesh = CuboidMesh(1, 1, 1, 1, 1, 1)
    f = scalar_from_coordinates(lambda r: r[0] + r[1] + r[2], mesh)
    assert allclose(f, np.array((1.5)))
