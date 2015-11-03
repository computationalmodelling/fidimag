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
                    np.array(((width / 2.0, height / 2.0, 0),
                              (width * 3.0 / 2.0, height / 2.0, 0))))


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
        (width / 2.0, height / 2.0, 0),
        (width * 3.0 / 2.0, height / 2.0, 0),
        (width, height * 5.0 / 4.0, 0),
        (width * 2, height * 5.0 / 4.0, 0))))


def test_neighbours_x():
    """

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    mesh = HexagonalMesh(1, 2, 1)
    # assert to_sets(mesh.neighbours) == [{1}, {0}]
    assert (mesh.neighbours == [[1, -1, -1, -1, -1, -1],
                                [-1, 0, -1, -1, -1, -1]]).all()


def test_neighbours_x_periodic():
    """

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, False))
    assert (mesh.neighbours == [[1, 1, -1, -1, -1, -1],
                                [0, 0, -1, -1, -1, -1]]).all()


def test_neighbours_x_periodic_all():
    """

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    print ""
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, True))
    assert (mesh.neighbours == [[1, 1, -1, -1, -1, -1],
                                [0, 0, -1, -1, -1, -1]]).all()


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
    # assert to_sets(mesh.neighbours) == [{1}, {0}]
    assert (mesh.neighbours == [[-1, -1, 1, -1, -1, -1],
                                [-1, -1, -1, 0, -1, -1]]).all()


def test_neighbours_y_square():
    """
    /\
   |  |
   | 1|
    \/\
    |  |    Hexagon size 1.
    | 0|    Cells 1 x 2.
     \/

    """
    mesh = HexagonalMesh(1, 1, 2, alignment='square')
    assert (mesh.neighbours == [[-1, -1, -1, -1, 1, -1],
                                [-1, -1, -1, -1, -1, 0]]).all()


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
    assert (mesh.neighbours == [[-1, -1, 1, 1, -1, -1],
                                [-1, -1, 0, 0, -1, -1]]).all()


# We need to CHECK the validity of the periodicity for this
# particular arrangement
# def test_neighbours_y_periodic_square():
#     """
#     /\
#    |  |
#    | 1|
#     \/\
#     |  |    Hexagon size 1.
#     | 0|    Cells 1 x 2.
#      \/
#
#     """
#     mesh = HexagonalMesh(1, 1, 2, periodicity=(False, True),
#                          alignment='square')
#     assert (mesh.neighbours == [[-1, -1, -1, -1, 1, 1],
#                                 [-1, -1, -1, -1, 0, 0]]).all()


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
    # assert to_sets(mesh.neighbours) == [{1, 3}, {0, 2, 3, 4},  # cells 0, 1
    #                                     {1, 4, 5}, {0, 1, 4, 6},  # cells 2, 3
    #                                     {1, 2, 3, 5, 6, 7}, {2, 4, 7, 8},  # cells 4, 5
    #                                     {3, 4, 7}, {4, 5, 6, 8},  # cells 6, 7
    #                                     {5, 7}]  # cell 8
    assert (mesh.neighbours == [[1, -1, 3, -1, -1, -1],  # cell 0
                                [2, 0, 4, -1, 3, -1],    # cell 1
                                [-1, 1, 5, -1, 4, -1],   # cell 2
                                [4, -1, 6, 0, -1, 1],    # cell 3
                                [5, 3, 7, 1, 6, 2],      # cell 4
                                [-1, 4, 8, 2, 7, -1],    # cell 5
                                [7, -1, -1, 3, -1, 4],   # ...
                                [8, 6, -1, 4, -1, 5],
                                [-1, 7, -1, 5, -1, -1]
                                ]
            ).all()


def test_neighbours_multiple_square():
    """

       /\ /\ /\
      |  |  |  |
      | 6| 7| 8|
     /\ /\ /\ /
    |  |  |  |
    | 3| 4| 5|
     \ /\ /\ /\
      |  |  |  |
      | 0| 1| 2|   Hexagon size 1.
       \/ \/ \/    Cells 3 x 3.

    """
    mesh = HexagonalMesh(1, 3, 3, alignment='square')
    print mesh.neighbours
    assert (mesh.neighbours == [[1, -1, 4, -1, 3, -1],   # cell 0
                                [2, 0, 5, -1, 4, -1],    # cell 1
                                [-1, 1, -1, -1, 5, -1],  # cell 2
                                [4, -1, 6, -1, -1, 0],   # ...
                                [5, 3, 7, 0, 6, 1],
                                [-1, 4, 8, 1, 7, 2],
                                [7, -1, -1, 3, -1, 4],
                                [8, 6, -1, 4, -1, 5],
                                [-1, 7, -1, 5, -1, -1]
                                ]
            ).all()

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
