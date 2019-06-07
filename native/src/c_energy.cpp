#include "c_energy.h"
#include "c_constants.h"

double Energy::compute_energy() {
    double sum = 0;
    for(int i = 0; i < n; i++) {
        sum += energy[i];
    }
    return sum * (dx * dy * dz * std::pow(unit_length, 3));
}


void Energy::setup(int nx, int ny, int nz, double dx, double dy, double dz,
                   double unit_length, double *spin, double *Ms, double *Ms_inv,
                   double *coordinates, double *ngbs, 
                   double *energy, double *field
                   ); {
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->unit_length = unit_length;

    // Arrays
    this->spin = spin;
    this->Ms = Ms;
    this->Ms_inv = Ms_inv;
    this->coordinates = coordinates;
    this->ngbs = ngbs;
    this ->energy->energy;
    this->field = field;

    set_up = true;
}

void ExchangeEnergy::compute_field(double t) {
    /* Compute the micromagnetic exchange field and energy using the
     * matrix of neighbouring spins and a second order approximation
     * for the derivative
     *
     * Ms_inv     :: Array with the (1 / Ms) values for every mesh node.
     *               The values are zero for points with Ms = 0 (no material)
     *
     * A          :: Exchange constant
     *
     * dx, dy, dz :: Mesh spacings in the corresponding directions
     *
     * n          :: Number of mesh nodes
     *
     * ngbs       :: The array of neighbouring spins, which has (6 * n)
     *               entries. Specifically, it contains the indexes of
     *               the neighbours of every mesh node, in the following order:
     *                      -x, +x, -y, +y, -z, +z
     *
     *               Thus, the array is like:
     *              | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  |
     *                i=0                           i=1                ...
     *
     *              where  0-y  is the index of the neighbour of the 0th spin,
     *              in the -y direction, for example. The index value for a
     *              neighbour where Ms = 0, is evaluated as -1. The array
     *              automatically gives periodic boundaries.
     *
     *              A basic example is a 3 x 3 two dimensional mesh with PBCs
     *              in the X and Y direction:
     *
     *                                  +-----------+
     *                                  | 6 | 7 | 8 |
     *                                  +-----------+
     *                                  | 3 | 4 | 5 |
     *                                  +-----------+
     *                                  | 0 | 1 | 2 |
     *                                  +-----------+
     *
     *              so, the first 6 entries (neighbours of the 0th mesh node)
     *              of the array would be: [ 2  1  6  3 -1 -1 ... ]
     *              (-1 since there is no material in +-z, and a '2' first,
     *              since it is the left neighbour which is the PBC in x, etc..)
     *
     *  For the exchange computation, the field is defined as:
     *          H_ex = (2 * A / (mu0 * Ms)) * nabla^2 (mx, my, mz)
     *
     *  Therefore, for the i-th mesh node (spin), we approximate the
     *  derivatives as:
     *          nabla^2 mx = (1 / dx^2) * ( spin[i-x] - 2 * spin[i] + spin[i+x] ) +
     *                       (1 / dy^2) * ( spin[i-y] - 2 * spin[i] + spin[i+y] ) +
     *                       (1 / dz^2) * ( spin[i-z] - 2 * spin[i] + spin[i+z] )
     *
     *  Where i-x is the neighbour in the -x direction. This is similar
     *  for my and mz.
     *  We can notice that the sum is the same if we do:
     *        ( spin[i-x] - spin[i] ) + ( spin[i+x] - spin[i] )
     *  so we can iterate through the neighbours and perform the sum with the
     *  corresponding coefficient 1 /dx, 1/dy or 1/dz
     *
     *  The *m array contains the spins as:
     *      [mx0, my0, mz0, mx1, my1, mz1, mx2, ...]
     *  so if we want the starting position of the magnetisation for the
     *  i-th spin, we only have to do (3 * i) for mx, (3 * i + 1) for my
     *  and (3 * i + 2) for mz
     *
     *
     *  IMPORTANT: The ex field usually has the structure:
     *                   2 * A / (mu0 Ms ) * (Second derivative of M)
     *       When discretising the derivative, it carries a "2" in the
     *       denominator which we "cancel" with the "2" in the prefactor,
     *       hence we do not put it explicitly in the calculations
     *
     *       So, when computing the energy: (-1/2) * mu * Ms * H_ex
     *       we only put the 0.5 factor and don't worry about the "2"s in the
     *       field
     *
     */

    /* Here we iterate through every mesh node */
	for (int i = 0; i < n; i++) {
        /* Define the coefficients */
        double ax = 2 * A[i] / (dx * dx);
        double ay = 2 * A[i] / (dy * dy);
        double az = 2 * A[i] / (dz * dz);

	    double fx = 0, fy = 0, fz = 0;
	    int idnm = 0;     // Index for the magnetisation matrix
	    int idn = 6 * i; // index for the neighbours

        /* Set a zero field for sites without magnetic material */
	    if (Ms_inv[i] == 0.0){
	        field[3 * i] = 0;
	        field[3 * i + 1] = 0;
	        field[3 * i + 2] = 0;
	        continue;
	    }

        /* Here we iterate through the neighbours */
        for (int j = 0; j < 6; j++) {
            /* Remember that index=-1 is for sites without material */
	        if (ngbs[idn + j] >= 0) {
                /* Magnetisation of the neighbouring spin since ngbs gives
                 * the neighbour's index */
	            idnm = 3 * ngbs[idn + j];

                /* Check that the magnetisation of the neighbouring spin
                 * is larger than zero */
                if (Ms_inv[ngbs[idn + j]] > 0){

                    /* Neighbours in the -x and +x directions
                     * giving: ( spin[i-x] - spin[i] ) + ( spin[i+x] - spin[i] )
                     * when ngbs[idn + j] > 0 for j = 0 and j=1
                     * If, for example, there is no
                     * neighbour at -x (j=0) in the 0th node (no PBCs),
                     * the second derivative would only be avaluated as:
                     *      (1 / dx * dx) * ( spin[i+x] - spin[i] )
                     * which, according to
                     * [M.J. Donahue and D.G. Porter; Physica B, 343, 177-183 (2004)]
                     * when performing the integration of the energy, we still
                     * have error of the order O(dx^2)
                     * This same applies for the other directions
                     */
                    if (j == 0 || j == 1) {
                        fx += ax * (spin[idnm]     - spin[3 * i]);
                        fy += ax * (spin[idnm + 1] - spin[3 * i + 1]);
                        fz += ax * (spin[idnm + 2] - spin[3 * i + 2]);
                    }
                    /* Neighbours in the -y and +y directions */
                    else if (j == 2 || j == 3) {
                        fx += ay * (spin[idnm]     - spin[3  * i]);
                        fy += ay * (spin[idnm + 1] - spin[3 * i + 1]);
                        fz += ay * (spin[idnm + 2] - spin[3 * i + 2]);
                    }
                    /* Neighbours in the -z and +z directions */
                    else if (j == 4 || j == 5) {
                        fx += az * (spin[idnm]     - spin[3 * i]);
                        fy += az * (spin[idnm + 1] - spin[3 * i + 1]);
                        fz += az * (spin[idnm + 2] - spin[3 * i + 2]);
                    }
                    else {
                        continue;
                    }
                }
            }
        }

        /* Energy as: (-mu0 * Ms / 2) * [ H_ex * m ]   */
        energy[i] = -0.5 * (fx * spin[3 * i] + fy * spin[3 * i + 1]
                            + fz * spin[3 * i + 2]);

        /* Update the field H_ex which has the same structure than *m */
        field[3 * i]     = fx * Ms_inv[i] * MU_0_INV;
        field[3 * i + 1] = fy * Ms_inv[i] * MU_0_INV;
        field[3 * i + 2] = fz * Ms_inv[i] * MU_0_INV;
    }
}
