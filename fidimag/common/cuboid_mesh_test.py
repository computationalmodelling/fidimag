import numpy as np
from cuboid_mesh import CuboidMesh


def allclose(a, b):
    """ close to machine precision """
    return np.allclose(a, b, rtol=1e-14, atol=1e-14)


def to_sets(arr):
    """ make sets out of the entries of arr for easier comparison """
    list_of_sets = [set(e) for e in list(arr)]
    for s in list_of_sets:
        if -1 in s:
            s.remove(-1)
    return list_of_sets


def test_coordinates_x():
    """
       .+---+---+
     .'  .'   .'|
    +---+---+'  |        Dimensions 2 x 1 x 1.
    |   |   |   |        Cells 2 x 1 x 1.
    |   |   |   +
    |   |   | .'
    +---+---+'

    """
    mesh = CuboidMesh(1, 1, 1, 2, 1, 1)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.5, 0.5), (1.5, 0.5, 0.5))))


def test_coordinates_y():
    """
       .+-------+
     .+-------+'|
    +-------+'| |        Dimensions 1 x 2 x 1.
    |       | | |        Cells 1 x 2 x 1.
    |       | | +
    |       | +'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 2, 1)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.5, 0.5), (0.5, 1.5, 0.5))))


def test_coordinates_z():
    """
       .+-------+
     .'       .'|
    +-------+'  +        Dimensions 1 x 1 x 2.
    |       | .'|        Cells 1 x 1 x 2.
    +-------+'  +
    |       | .'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 1, 2)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.5, 0.5), (0.5, 0.5, 1.5))))


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
    assert to_sets(mesh.neighbours) == [{1}, {0}]


def test_neighbours_x_periodic():
    """
       .+---+---+
     .'  .'   .'|
    +---+---+'  |
    |   |   |   |
    | 0 | 1 |   +
    |   |   | .'
    +---+---+'

    """
    mesh = CuboidMesh(1, 1, 1, 2, 1, 1, periodicity=(True, False, False))
    # over under behind in-front right left
    assert allclose(mesh.neighbours, np.array(((1, 1, -1, -1, -1, -1), (0, 0, -1, -1, -1, -1))))


def test_neighbours_x_periodic_all():
    """
       .+---+---+
     .'  .'   .'|
    +---+---+'  |
    |   |   |   |
    | 0 | 1 |   +
    |   |   | .'
    +---+---+'

    """
    mesh = CuboidMesh(1, 1, 1, 2, 1, 1, periodicity=(True, True, True))
    assert allclose(mesh.neighbours, np.array(((1, 1, -1, -1, -1, -1), (0, 0, -1, -1, -1, -1))))


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
    assert to_sets(mesh.neighbours) == [{1}, {0}]


def test_neighbours_y_periodic():
    """
       .+-------+
     .+-------+'|
    +-------+'| |
    |       | |1|
    |       |0| +
    |       | +'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 2, 1, periodicity=(False, True, False))
    assert allclose(mesh.neighbours, np.array(((-1, -1, 1, 1, -1, -1), (-1, -1, 0, 0, -1, -1))))


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
    assert to_sets(mesh.neighbours) == [{1}, {0}]


def test_neighbours_z_periodic():
    """
       .+---0---+
     .'       .'|
    +-------+'  +
    |   1   | .'|
    +-------+'  +
    |   0   | .'
    +-------+'

    """
    mesh = CuboidMesh(1, 1, 1, 1, 1, 2, periodicity=(False, False, True))
    assert allclose(mesh.neighbours, np.array(((-1, -1, -1, -1, 1, 1), (-1, -1, -1, -1, 0, 0))))


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
    assert to_sets(mesh.neighbours) == [{1, 2, 4}, {0, 3, 5},  # for cells 0, 1
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
