#include "m_clib.h"
#include "c_vectormath.h"


// Compute: S \cdot (S_i \times S_j)
inline double volume(double S[3], double Si[3], double Sj[3]) {
	double tx = S[0] * (-Si[2] * Sj[1] + Si[1] * Sj[2]);
	double ty = S[1] * (Si[2] * Sj[0] - Si[0] * Sj[2]);
	double tz = S[2] * (-Si[1] * Sj[0] + Si[0] * Sj[1]);
	return tx + ty + tz;
}

double skyrmion_number(double * spin, double * charge,
                       int nx, int ny, int nz, int * ngbs) {

    /* Compute the skyrmion number Q, defined as:
     *                      _
     *           1         /       dM     dM
     *   Q   =  ---  *    /   M .  --  X  --   dx dy
     *          4 PI   _ /         dx     dy
     *
     *
     * for a two dimensional layer (it can be a slice of a
     * bulk system)
     *
     * Therefore, a finite differences discretisation of the
     * continuum magnetisation field, using central
     * differences, leads to:
     *
     * Q =   M_i \cdot ( M_{i+1} \times M_{j+1} )
     *       + M_i \cdot ( M_{i-1} \times M_{j-1} )
     *       - M_i \cdot ( M_{i-1} \times M_{j+1} )
     *       - M_i \cdot ( M_{i+1} \times M_{j-1} )
     *
     *
     * INPUTS:
     *
     * The *spin array is the vector field for a two dimensional
     * lattice with dimensions nx * ny
     * (we can take a slice of a bulk from Python and pass it here,
     *  remember to do the ame for the neighbours matrix)
     * The array follows the order:
     *   [Sx0 Sy0 Sz0 Sx1 Sy1 Sz1 ... ]
     *
     * Charge is a scalar field array used to store the spin chirality /
     * skyrmion number density (Q value per lattice site)
     *
     * *ngbs is the array with the neighbours information for every
     * lattice site. The array is like:
     *      [ 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  ]
     *        i=0                           i=1                ...
     *
     * where  0-y  is the index of the neighbour of the 0th spin,
     * in the -y direction, for example
     */

	int i, j;
	int index, id;

	double sum = 0;

	double S[3];

	int nxy = nx * ny;

    /* Store the spin directions of the nearest neighbours
     * in the order: [-x +x -y +y]
     */
    double S_nn[12];

        for(i=0;i<12;i++){
          S_nn[i] = 0; //we have to set S_nn to zeros manually
        }

	for (i = 0; i < nxy; i++) {
        index = 3 * i;

        /* The starting index of the nearest neighbours for the
         * i-th spin */
        int id_nn = 6 * i;

        S[0] = spin[index];
        S[1] = spin[index + 1];
        S[2] = spin[index + 2];

        for (j = 0; j < 4; j++) {
            if (ngbs[id_nn + j] > 0) {
                id = 3 * ngbs[id_nn + j];
                S_nn[3 * j    ] = spin[id];
                S_nn[3 * j + 1] = spin[id + 1];
                S_nn[3 * j + 2] = spin[id + 2];
            }
        }

        charge[i] =   volume(S, &S_nn[3], &S_nn[9])
                    + volume(S, &S_nn[0], &S_nn[6])
                    - volume(S, &S_nn[0], &S_nn[9])
                    - volume(S, &S_nn[3], &S_nn[6]);

        charge[i] /= (16 * M_PIl);

        sum += charge[i];
    }

	return sum;

}
