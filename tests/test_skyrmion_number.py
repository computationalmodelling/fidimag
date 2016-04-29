from __future__ import print_function
from fidimag.atomistic import Sim
from fidimag.common import CuboidMesh
from fidimag.atomistic import DMI
from fidimag.atomistic import UniformExchange
from fidimag.atomistic import Zeeman

from fidimag.micro import DMI as microDMI
from fidimag.micro import UniaxialAnisotropy as microUniaxialAnisotropy
from fidimag.micro import UniformExchange as microUniformExchange
from fidimag.micro import Zeeman as microZeeman
from fidimag.micro import Sim as microSim


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

    This test generate a skyrmion pointing down with unrealistic
    paremeters.

    """

    mesh = CuboidMesh(nx=120, ny=120, nz=1,
                      periodicity=(True, True, False))

    sim = Sim(mesh, name='skx_num')
    sim.set_tols(rtol=1e-6, atol=1e-6)
    sim.alpha = 1.0
    sim.gamma = 1.0
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
    assert skn > -1 and skn < -0.99


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
    sim.set_tols(rtol=1e-6, atol=1e-6)
    sim.alpha = 1.0
    sim.gamma = 1.0
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
    test_skx_num_micromagnetic()
