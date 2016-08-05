import unittest


def regression_helper_init_scalar():
    import fidimag
    import numpy as np
    # Create simulation
    mesh = fidimag.common.CuboidMesh(nx=10, ny=2, nz=1,
                                     dx=1, dy=1, dz=1)
    sim = fidimag.atomistic.Sim(mesh)
    # Set the magnetisation and save it to mu_s.npy
    sim.mu_s = 3
    np.save('mu_s.npy', sim.mu_s)

    # Create another simulation with less elements in the x direction
    mesh2 = fidimag.common.CuboidMesh(nx=5, ny=2, nz=1,
                                      dx=1, dy=1, dz=1)
    sim2 = fidimag.atomistic.Sim(mesh2)
    # Set the magnetisation loading the other simulation file
    sim2.mu_s = np.load('mu_s.npy')


class regression_test(unittest.TestCase):

    def test_init_scalar(self):
        self.assertRaises(ValueError, regression_helper_init_scalar)
