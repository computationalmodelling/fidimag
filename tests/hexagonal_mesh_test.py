from __future__ import print_function
import pytest
import numpy as np
from math import sqrt
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh


def allclose(a, b):
    """ close to machine precision """
    return np.allclose(a, b, rtol=1e-14, atol=1e-14)


def to_sets(arr):
    """ make sets out of the entries of arr for easier comparison """
    return [set(e) for e in arr]


def test_coordinates_x():
    r"""

     /\ /\
    |  |  |  Hexagon size 1.
    |  |  |  Cells 2 x 1.
     \/ \/


    """
    size = 1
    mesh = HexagonalMesh(size, 2, 1)
    width = size * 2.0
    height = 2.0 * width / sqrt(3)
    assert allclose(mesh.coordinates,
                    np.array(((width / 2.0, height / 2.0, 0),
                              (width * 3.0 / 2.0, height / 2.0, 0))))


def test_coordinates_y():
    r"""
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
    width = size * 2.0
    # height refers to the y-distance between two hexagons
    # centres in consecutive rows
    height = sqrt(3) * size
    # This is the total hexagon height that we use as the base
    hex_height = 2.0 * width / sqrt(3)
    mesh = HexagonalMesh(size, 2, 2)
    assert allclose(mesh.coordinates, np.array((
        (width / 2.0, hex_height / 2.0, 0),
        (width * 3.0 / 2.0, hex_height / 2.0, 0),
        (width, height + hex_height / 2.0, 0),
        (width * 2, height + hex_height / 2.0, 0))))


def test_neighbours_x():
    r"""

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
    r"""

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, False))
    assert (mesh.neighbours == [[1, 1, -1, -1, -1, -1],
                                [0, 0, -1, -1, -1, -1]]).all()


def test_neighbours_x_periodic_all():
    r"""

     /\ /\
    |  |  |  Hexagon size 1.
    | 0| 1|  Cells 2 x 1.
     \/ \/

    """
    print("")
    mesh = HexagonalMesh(1, 2, 1, periodicity=(True, True))
    assert (mesh.neighbours == [[1, 1, -1, -1, -1, -1],
                                [0, 0, -1, -1, -1, -1]]).all()


def test_neighbours_y():
    r"""
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
    r"""
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
    r"""
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


def test_nearest_neighbours_multiple():
    r"""
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
    print(mesh.neighbours)
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


def test_nearest_neighbours_multiple_square():
    r"""

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
    print(mesh.neighbours)
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


# -----------------------------------------------------------------------------

def test_neighbours_9shells_square_odd_row():
    r"""

    Testing neighbour indexes manually, for a 11x11 hexagonal mesh with square
    alignment. In this test we check the neighbours of the 60th lattice site,
    which is in an even row (we could check an odd row as well in the future)

       /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
      |  |  |  |  |..|..|  |  |  |  |  |
      |  |  |  |  |14|15|  |  |  |  |  |
     /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /
    |  |  |  |;;|~~|``|~~|;;|  |  |  |
    |  |  |  |02|03|04|05|06|  |  |  |
     \ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
      |  |..|~~|##|@@|@@|##|~~|..|  |  |
      |  |89|90|91|92|93|94|95|96|  |  |
     /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /
    |  |..|``|@@|**|--|**|@@|``|..|  |
    |  |78|79|80|81|82|83|84|85|86|  |
     \ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
      |  |~~|@@|--|xx|xx|--|@@|~~|  |  |          OO i-th site
      |  |67|68|69|70|71|72|73|74|  |  |          xx 1st neighbours (nearest)
     /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /           -- 2nd neighbours (next nearest)
    |  |;;|##|**|xx|  |xx|**|##|;;|  |            ** 3rd neighbours
    |  |56|57|58|59|60|61|62|63|64|  |            @@ 4th neighbours
     \ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\           ## 5th neighnbours
      |  |~~|@@|--|xx|xx|--|@@|~~|  |  |          `` 6th neighbours
      |  |45|46|47|48|49|50|51|52|  |  |          ~~ 7th neighbours
     /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /           ;; 8th neighbours
    |  |..|``|@@|**|--|**|@@|``|..|  |
    |  |34|35|36|37|38|39|40|41|42|  |
     \ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
      |  |..|~~|##|@@|@@|##|~~|..|  |  |
      |22|23|24|25|26|27|28|29|30|  |  |
     /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /
    |  |  |  |;;|~~|``|~~|;;|  |  |  |
    |11|  |13|14|15|16|17|18|  |  |  |
     \ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
      |  |  |  |  |..|..|  |  |  |  |  |
      | 0|  | 2|  | 4| 5|  |  |  |  |  |         Hexagon size 1.
       \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/          Cells 9 x 6.


    """
    mesh = HexagonalMesh(1, 11, 11, alignment='square', shells=9)
    print(mesh.neighbours)
    assert (mesh.neighbours[60] == [61, 59, 71, 48, 70, 49,    # 1st shell
                                    72, 47, 82, 38, 69, 50,    # 2nd shell
                                    62, 58, 83, 37, 81, 39,    # 3rd shell
                                    73, 46, 84, 36, 93, 26,    # 4th shell
                                    92, 27, 80, 40, 68, 51,
                                    63, 57, 94, 25, 91, 28,    # 5th shell
                                    85, 35, 104, 16, 79, 41,   # 6th shell
                                    74, 45, 95, 24, 105, 15,   # 7th shell
                                    103, 17, 90, 29, 67, 52,
                                    64, 56, 106, 14, 102, 18,  # 8th shell
                                    86, 34, 96, 23, 115, 4,    # 9th shell
                                    114, 5, 89, 30, 78, 42
                                    ]
            ).all()

    # Second neighbours of the 0th position:
    assert (mesh.neighbours[0][6:12] == [13, -1, 22, -1, -1, -1]  # 2nd shell
            ).all()


def test_neighbours_9shells_square_even_row():
    """
    Just as test_neighbours_9shells_square_odd_row but checking a site in an
    even row in a 12x12 hexagonal lattice
    """
    mesh = HexagonalMesh(1, 12, 12, alignment='square', shells=9)
    assert (mesh.neighbours[77] == [78, 76, 90, 65, 89, 66,    # 1st shell
                                    91, 64, 101, 53, 88, 67,   # 2nd shell
                                    79, 75, 102, 52, 100, 54,  # 3rd shell
                                    92, 63, 103, 51, 114, 41,  # 4th shell
                                    113, 42, 99, 55, 87, 68,
                                    80, 74, 115, 40, 112, 43,  # 5th shell
                                    104, 50, 125, 29, 98, 56,  # 6th shell
                                    93, 62, 116, 39, 126, 28,  # 7th shell
                                    124, 30, 111, 44, 86, 69,
                                    81, 73, 127, 27, 123, 31,  # 8th shell
                                    105, 49, 117, 38, 138, 17, # 9th shell
                                    137, 18, 110, 45, 97, 57
                                    ]
            ).all()


def test_iterate_over_cells():
    mesh = HexagonalMesh(1, 2, 2)
    for c_i in mesh.cells():
        print("This is cell #{}, I have neighbours {}.".format(c_i, mesh.neighbours[c_i]))


def test_iterate_over_cells_and_neighbours():
    mesh = HexagonalMesh(1, 2, 2)
    for c_i in mesh.cells():
        print("I am cell #{}.".format(c_i))
        for c_j in mesh.neighbours[c_i]:
            print("\tAnd I am its neighbour, cell #{}!".format(c_j))


@pytest.mark.xfail(reason="Skipping because this is not supported")
def test_hexagonal_mesh_creation_periodic_x():
    mesh = HexagonalMesh(1, 2, 2, alignment='square', periodicity=(True, False, False))


@pytest.mark.xfail(reason="Skipping because this is not supported")
def test_hexagonal_mesh_creation_periodic_square_y():
    mesh = HexagonalMesh(1, 2, 2, alignment='square', periodicity=(False, True, False))
