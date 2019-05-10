import numpy as np
import fidimag.extensions.c_clib as clib


def test_sst_field_1d():
    r"""
    This is a direct test of the STT C library
    We create a 1-D  4 spins system along the x direction with the
    following components:

               s_x   s_y   s_z       spin
            [[  0.    4.    8. ]    i = 0
             [  0.8   5.    9. ]    i = 1
             [  2.    6.   10. ]    i = 2
             [  2.1   7.   11. ]]   i = 3

    Most of the STT parameters are set to zero, so the field is reduced
    to the calculation of:

    H_STT =  j \cdot \nabla S = ( jx * d S_x / dx + jy * d S_x / dy ,
                                  jx * d S_y / dx + jy * d S_y / dy ,
                                  0)

    Therefore, the x component of the field for the 1th spin, for example,
    would be:

            H_STTx (i=1) = jx * d S_x (i=1) / dx + 0
                         = [ S_x (i=2)  -  S_x (i=0) ]  / 2 * dx
                         = ( 2 - 0 ) / 2
                         = 1

    since jx = dx = 1 and the derivatives along y are zero (1D spin chain)

    The first test is the spin chain without PBC and the second
    test is using PBCs, so the neighbours matrix changes

    ** Without PBC:
    At the extremes, the derivative has a different discretisation
    (using only the right or left neighbour and centered at the spin)
    which, for the 0th spin, for instance, is
            H_STTx (i=0) = jx * d S_x (i=0) / dx + 0
                         = [ S_x (i=1)  -  S_x (i=0) ]  / dx
                         = ( 0.8 - 0 )
                         = 0.8

    For further details, see: fidimag/atomistic/lib/stt.c

    With PBC: The derivatives at the extremes are performed with the
    periodic neighbour, e.g. for the 0th spin, the NN at -x is the 3th spin

    """

    # The system is 1d with 4 spins.
    nx, ny, nz = 4, 1, 1
    n = nx * ny * nz
    dx, dy, dz = 1, 1, 1

    # NO PBCs -----------------------------------------------------------------
    # We can manually construct the neighbours matrix,
    # the neighbours order is [-x, +x, -y, +y, -z, +z]
    # Since there are no PBCs, for the 0th spin, i.e. ngbs[0],
    # the left (-x) NN is -1 and for the 3th spin, the +x NN is -1
    # Also, since it is a 1D chain, there are no NNs in y or z
    ngbs = [[-1, 1, -1, -1, -1, -1],
            [0, 2, -1, -1, -1, -1],
            [1, 3, -1, -1, -1, -1],
            [2, -1, -1, -1, -1, -1],
            ]
    ngbs = np.array(ngbs, dtype=np.int32)

    nxyz = nx * ny * nz * 3

    # The spin values to get the matrix:
    #
    #         [[  0.   4.    8.  ]
    #          [  1    5.    9.  ]
    #          [  2.   6.    10. ]
    #          [  3    7.    11. ]]
    spin = np.array([1.0 * i for i in range(nxyz)]).reshape((-1, 3), order='F').flatten()

    field = np.zeros(nxyz)
    jx = np.zeros(nxyz)
    jy = np.zeros(nxyz)
    jz = np.zeros(nxyz)
    jx[:] = 1
    jy[:] = 2

    # Modify the x components of the spins to get the values
    # specified in the problem description
    spin[3 * 1] = 0.8
    spin[3 * 2] = 2
    spin[3 * 3] = 2.1

    # print spin.reshape(-1, 3)

    clib.compute_stt_field(spin,
                           field,
                           jx, jy, jz,
                           dx, dy, dz, 
                           ngbs,
                           n
                           )
    # print field

    assert field[3 * 0] == 0.8
    assert field[3 * 1] == 1
    assert field[3 * 2] == 0.65
    assert abs(field[3 * 3] - 0.1) < 1e-16

    # PBCs --------------------------------------------------------------------
    # Now we make PBCs in the x direction so the -x NN for
    # the 0th spin is the 3rd spin, and the opposite way for the 3rd spin
    ngbs = [[3, 1, -1, -1, -1, -1],
            [0, 2, -1, -1, -1, -1],
            [1, 3, -1, -1, -1, -1],
            [2, 0, -1, -1, -1, -1],
            ]
    ngbs = np.array(ngbs, dtype=np.int32)

    clib.compute_stt_field(spin,
                           field,
                           jx, jy, jz,
                           dx, dy, dz,
                           ngbs,
                           n
                           )
    # print field
    assert field[3 * 0] == -0.65
    assert field[3 * 3] == -1


if __name__ == '__main__':
    test_sst_field_1d()
