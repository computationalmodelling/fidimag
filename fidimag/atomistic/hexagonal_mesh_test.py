import numpy as np
from math import sqrt
from hexagonal_mesh import HexagonalMesh


def allclose(a, b):
    """ close to machine precision """
    return np.allclose(a, b, rtol=1e-14, atol=1e-14)


def to_sets(arr):
    """ make sets out of the entries of arr for easier comparison """
    return [set(e) for e in arr]


def test_coordinates_x():
    """

     /\ /\
    |  |  |  Hexagon size 1.
    |  |  |  Cells 2 x 1.
     \/ \/


    """
    size = 1
    mesh = HexagonalMesh(size, 2, 1)
    height = size * 2.0
    width = sqrt(3) / 2.0 * height
    assert allclose(mesh.coordinates,
                    np.array(((width / 2.0, height / 2.0),
                              (width * 3.0 / 2.0, height / 2.0))))


def test_coordinates_y():
    """
       /\ /\
      |  |  |
      | 2| 3|
     /\ /\ /
    |  |  |    Hexagon size 1.
    | 0| 1|    Cells 2 x 2.
     \/ \/

    """
    size = 1
    mesh = HexagonalMesh(size, 2, 1)
    height = size * 2.0
    width = sqrt(3) / 2.0 * height
    mesh = HexagonalMesh(size, 2, 2)
    assert allclose(mesh.coordinates, np.array((
        (width / 2.0, height / 2.0),
        (width * 3.0 / 2.0, height / 2.0),
        (width, height * 5.0 / 4.0),
        (width * 2, height * 5.0 / 4.0))))


def test_neighbours_x():
    """

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    mesh = HexagonalMesh(1, 2, 1)
    assert to_sets(mesh.neighbours) == [{1}, {0}]


def test_neighbours_x_periodic():
    """
 
     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, False))
    assert mesh.neighbours == [[1, 1], [0, 0]]


def test_neighbours_x_periodic_all():
    """
 
     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    print ""
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, True))
    assert mesh.neighbours == [[1, 1], [0, 0]]


def test_neighbours_y():
    """
       /\ 
      |  |
      | 1|
     /\ /
    |  |    Hexagon size 1.
    | 0|    Cells 1 x 2.
     \/ 

    """
    mesh = HexagonalMesh(1, 1, 2)
    assert to_sets(mesh.neighbours) == [{1}, {0}]


def test_neighbours_y_periodic():
    """
       /\ 
      |  |
      | 1|
     /\ /
    |  |    Hexagon size 1.
    | 0|    Cells 1 x 2.
     \/ 

    """
    mesh = HexagonalMesh(1, 1, 2, periodicity=(False, True))
    assert mesh.neighbours == [[1, 1], [0, 0]]


def test_neighbours_multiple():
    """
         /\ /\ /\
        |  |  |  |
        | 6| 7| 8|
       /\ /\ /\ /
      |  |  |  |
      | 3| 4| 5|
     /\ /\ /\ /
    |  |  |  |    Hexagon size 1.
    | 0| 1| 2|    Cells 3 x 3.
     \/ \/ \/


    """
    mesh = HexagonalMesh(1, 3, 3)
    print mesh.neighbours
    assert to_sets(mesh.neighbours) == [{1, 3}, {0, 2, 3, 4},  # cells 0, 1
                                        {1, 4, 5}, {0, 1, 4, 6},  # cells 2, 3
                                        {1, 2, 3, 5, 6, 7}, {2, 4, 7, 8},  # cells 4, 5
                                        {3, 4, 7}, {4, 5, 6, 8},  # cells 6, 7
                                        {5, 7}]  # cell 8


def test_iterate_over_cells():
    mesh = HexagonalMesh(1, 2, 2)
    for c_i in mesh.cells():
        print "This is cell #{}, I have neighbours {}.".format(c_i, mesh.neighbours[c_i])


def test_iterate_over_cells_and_neighbours():
    mesh = HexagonalMesh(1, 2, 2)
    for c_i in mesh.cells():
        print "I am cell #{}.".format(c_i)
        for c_j in mesh.neighbours[c_i]:
            print "\tAnd I am its neighbour, cell #{}!".format(c_j)
