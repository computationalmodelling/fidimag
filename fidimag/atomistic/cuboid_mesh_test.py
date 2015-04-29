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
       .+---+---+
     .+-------+'|
    +---+---+'| |        Dimensions 1 x 1 x 1.
    |       | | |        Cells 1 x 2 x 1.
    |       | | +   
    |       | +'  
    +---+---+'   

    """
    mesh = CuboidMesh(1, 1, 1, 1, 2, 1)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.25, 0.5), (0.5, 0.75, 0.5))))


def test_coordinates_z():
    """
       .+---+---+
     .'       .'|
    +---+---+'  +        Dimensions 1 x 1 x 1.
    |       | .'|        Cells 1 x 1 x 2.
    +---+---+'  +   
    |       | .'  
    +---+---+'   

    """
    mesh = CuboidMesh(1, 1, 1, 1, 1, 2)
    assert allclose(mesh.coordinates,
                    np.array(((0.5, 0.5, 0.25), (0.5, 0.5, 0.75))))
