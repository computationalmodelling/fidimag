from __future__ import print_function
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import DMI
from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Zeeman
from fidimag.atomistic import Anisotropy

from fidimag.micro import DMI as microDMI
from fidimag.micro import UniaxialAnisotropy as microUniaxialAnisotropy
from fidimag.micro import UniformExchange as microUniformExchange
from fidimag.micro import Zeeman as microZeeman
from fidimag.micro import Sim as microSim

import fidimag.common.constant as const

import numpy as np


def init_m(pos, x0, y0, r):
    x = pos[0]
    y = pos[1]

    # x0, y0, r = 60, 60, 20

    m1 = (0.05, 0.01, -1)
    m2 = (0, 0, 1)

    if (x - x0) ** 2 + (y - y0) ** 2 < r ** 2:
        return m1
    else:
        return m2


def init_m_multiple_sks(pos, r, sk_pos):
    # Generate singularities with radius r at the positions in the
    # sk_pos list, in order to generate multiple skyrmions after
    # relaxation. Example: sk_pos=[(1, 1), (2, 2)]
    #                               x  y
    x, y = pos[0], pos[1]

    for coord in sk_pos:
        if (x - coord[0]) ** 2 + (y - coord[1]) ** 2 < r ** 2:
            return (0, 0, 1)

    return (0, 0, -1)


def test_skx_num_atomistic():
    """
    Test the *finite spin chirality* or skyrmion number for
    a discrete spins simulation in a two dimensional lattice

    The expression is (PRL 108, 017601 (2012)) :

    Q =     S_i \dot ( S_{i+1}  X  S_{j+1} )
         +  S_i \dot ( S_{i-1}  X  S_{j-1} )

    which measures the chirality taking two triangles of spins
    per lattice site i:
        S_{i} , S_{i + x} , S_{i + y}    and
        S_{i} , S_{i - x} , S_{i - y}

    The area of the two triangles cover a unit cell, thus the sum
    cover the whole area of the atomic lattice

    We also test the Berg and Luscher definition for a topological
    charge (see the hexagonal mesh test for details) in a
    square lattice.

    This test generate a skyrmion pointing down with unrealistic
    paremeters.

    """

    mesh = CuboidMesh(nx=120, ny=120, nz=1,
                      periodicity=(True, True, False))

    sim = Sim(mesh, name='skx_num')
    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 1.0
    sim.driver.gamma = 1.0
    sim.mu_s = 1.0

    sim.set_m(lambda pos: init_m(pos, 60, 60, 20))

    sim.do_precession = False

    J = 1.0
    exch = UniformExchange(J)
    sim.add(exch)

    D = 0.09
    dmi = DMI(D)
    sim.add(dmi)

    zeeman = Zeeman([0, 0, 5e-3])
    sim.add(zeeman)

    sim.relax(dt=2.0, stopping_dmdt=1e-2, max_steps=1000,
              save_m_steps=None, save_vtk_steps=None)

    skn = sim.skyrmion_number()
    print('skx_number', skn)

    skn_BL = sim.skyrmion_number(method='BergLuscher')
    print('skx_number_BergLuscher', skn_BL)

    # Test the finite chirality method
    assert skn > -1 and skn < -0.99

    # Test the Berg-Luscher method
    assert np.abs(skn_BL - (-1)) < 1e-4 and np.sign(skn_BL) < 0


def test_skx_num_atomistic_hexagonal():
    """

    Test the topological charge or skyrmion number for a discrete spins
    simulation in a two dimensional hexagonal lattice, using Berg and Luscher
    definition in [Nucl Phys B 190, 412 (1981)] and simplified in [PRB 93,
    174403 (2016)], which maps a triangulated lattice (using triangles of
    neighbouring spins) area into a unit sphere.

    The areas of two triangles per lattice site cover a unit cell, thus the sum
    cover the whole area of the atomic lattice

    This test generates a skyrmion pointing down and two skyrmions pointing up
    in a PdFe sample using magnetic parameters from: PRL 114, 177203 (2015)

    """

    mesh = HexagonalMesh(0.2715, 41, 41, periodicity=(True, True))

    sim = Sim(mesh, name='skx_number_hexagonal')
    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 1.0
    sim.driver.gamma = 1.0
    sim.mu_s = 3 * const.mu_B

    sim.set_m(lambda pos: init_m(pos, 16.1, 10, 2))

    sim.driver.do_precession = False

    J = 5.881 * const.meV
    exch = UniformExchange(J)
    sim.add(exch)

    D = 1.557 * const.meV
    dmi = DMI(D, dmi_type='interfacial')
    sim.add(dmi)

    sim.add(Anisotropy(0.406 * const.meV, axis=[0, 0, 1]))

    zeeman = Zeeman([0, 0, 2.5])
    sim.add(zeeman)

    sim.relax(dt=1e-13, stopping_dmdt=1e-2, max_steps=2000,
              save_m_steps=None, save_vtk_steps=100)

    skn_single = sim.skyrmion_number(method='BergLuscher')
    print('skx_number_hexagonal', skn_single)

    # Now we generate two skyrmions pointing up
    sim.driver.reset_integrator()
    sim.set_m(lambda pos: init_m_multiple_sks(pos, 1,
                                              sk_pos=[(9, 6), (18, 12)]
                                              )
              )
    sim.get_interaction('Zeeman').update_field([0, 0, -2.5])
    sim.relax(dt=1e-13, stopping_dmdt=1e-2, max_steps=2000,
              save_m_steps=None, save_vtk_steps=None)

    skn_two = sim.skyrmion_number(method='BergLuscher')
    print('skx_number_hexagonal_two', skn_two)

    # Check that we get a right sk number
    assert np.abs(skn_single - (-1)) < 1e-4 and np.sign(skn_single) < 0
    assert np.abs(skn_two - (2)) < 1e-4 and np.sign(skn_two) > 0


def test_skx_num_micromagnetic():
    """

    Compute the skyrmion number from a micromagnetic simulation,
    using the continuum topological charge expression for a
    2D layer:
                          _
               1         /       dM     dM
       Q   =  ---  *    /   M .  --  X  --   dx dy
              4 PI   _ /         dx     dy

    The skyrmion is generated in an Iron based sample, whose
    magnetic parameters are in: PRL 114, 177203 (2015)

    """

    # A 6 nm by 6 nm square sample with enough resolution for a
    # skyrmion with Q ~= -1
    mesh = CuboidMesh(nx=60, ny=60, nz=1, dx=0.1, dy=0.1, dz=0.4,
                      unit_length=1e-9,
                      periodicity=(True, True, False))

    sim = microSim(mesh, name='skx_num_micro')
    sim.driver.set_tols(rtol=1e-6, atol=1e-6)
    sim.driver.alpha = 1.0
    sim.driver.gamma = 1.0
    sim.Ms = 1.1e6

    sim.set_m(lambda pos: init_m(pos, x0=3, y0=3, r=1))

    sim.do_precession = False

    A = 2e-12
    exch = microUniformExchange(A)
    sim.add(exch)

    D = 3.9e-3
    dmi = microDMI(D, dmi_type='interfacial')
    sim.add(dmi)

    zeeman = microZeeman([0, 0, 2])
    sim.add(zeeman)

    Ku = 2.5e6
    sim.add(microUniaxialAnisotropy(Ku,
                                    (0, 0, 1), name='Ku'))

    sim.relax(stopping_dmdt=1e-2, max_steps=1000,
              save_m_steps=None,
              save_vtk_steps=None)

    sim.save_vtk()

    skn = sim.skyrmion_number()

    print('skx_number', skn)
    print(sim._skx_number)
    assert skn > -1 and skn < -0.99


if __name__ == '__main__':

    test_skx_num_atomistic()
    test_skx_num_atomistic_hexagonal()
    test_skx_num_micromagnetic()
