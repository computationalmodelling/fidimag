#include "a_clib.h"
#include "math.h"
#include "stdlib.h"

void demag_full(double * spin, double * field, double * energy, double * coords,
                double * mu_s, double * mu_s_scale, int n) {

    /* Full calculation of the Demagnetising field for atomistic systems
     * The main idea is to iterate through every lattice site, and sum
     * the dipolar contribution of the whole system. This is not the most
     * effective aproach since it needs N^2 calculations, but it is reliable
     * to compare against other approximations or techniques
     *
     * Thus, for the i-th spin, the field is calculated as:
     *
     *                               ^      ^         ^        ^
     *   ->      mu0 mu_s    __    3 r_ij ( m_j \cdot r_ij ) - m_j
     *   H_i =   --------   \	   -------------------------------
     *             4 pi     /__              r_ij ^ 3
     *
     *                    i != j
     *
     * where the numerator has unit vectors.
     * r_ij is a vector from r_i to r_j , i.e.  r_ij = r_j - r_i
     *
     * The prefactor for every lattice site, is stored in the mu_s_scale
     * array, since mu_s will depend on the site (different or no material)
     *
     * The spin and field are vector fields, thus they have 3 * n entries,
     * where n is the number of lattice sites of the system (including the
     * ones with mu_s = 0)
     * The coords array also has 3 * n entries since it has the coordinates
     * for every point, in the order:   x0 y0 z0 x1 y1 z1 x2 ...
     *
     * mu_s has the magnetic moments magnitudes and it is used for computing
     * the energy density of the i-th spin as:
     *
     *             mu_s    __   ^         ->
     *   E_i =  -   --    \     m_i \cdot H_i
     *               2    /__
     *
     *                  i=x,y,z
     *
     * the 1/2 factor is for the double counting.
     *
     */


    /* rij for the distance vector and rij_n for the normalised
     * version of the vector. rij_mag is the magnitude */
    double rij[3];
    double rij_n[3];

    /* we start iterating through every lattice site */
	#pragma omp parallel for private(rij, rij_n)
	for (int i = 0; i < n; i++) {
        double rij_mag;
        /* This is for the dot product of m and rij */
        double mrij = 0;

        /* Reset the field values */
        for (int k = 0; k < 3; k++) field[3 * i + k] = 0;

        /* Now we iterate through every spin of the system excluding
         * the point where we are right now */
        for (int j = 0; j < n; j++) {

            /* Avoid sites with mu_s = 0 */
            if(j != i && mu_s_scale[j] != 0.){

                /* Compute the distance from r_i to r_j and normalise */
                for(int k = 0; k < 3; k++) {
                    rij[k] = coords[3 * j + k] - coords[3 * i + k];
                }

                rij_mag = sqrt(rij[0] * rij[0] +
                               rij[1] * rij[1] +
                               rij[2] * rij[2]);

                for(int k = 0; k < 3; k++) rij_n[k] = rij[k] / rij_mag;

                /* dot product of m_j and r_ij (normalised)
                 * Remember that m has the structure: mx0 my0 mz0 mx1 my1 ...
                 */
                mrij = spin[3 * j] * rij_n[0] + spin[3 * j + 1] * rij_n[1]
                       + spin[3 * j + 2] * rij_n[2] ;

                /* Now add the contribution of the j-th spin */
                for(int k = 0; k < 3; k++){
                    field[3 * i + k] += (3 * rij_n[k] * mrij - spin[3 * j + k])
                                         / (rij_mag * rij_mag * rij_mag) ;
                }
            }
        }

        /* Now we scale the total contribution to the i-th spin */
        for(int k = 0; k < 3; k++) field[3 * i + k] *= mu_s_scale[i];

        /* And compute the energy avoiding double counting */
        energy[i] = -0.5 * mu_s[i] * (field[3 * i]     * spin[3 * i]     +
                                      field[3 * i + 1] * spin[3 * i + 1] +
                                      field[3 * i + 2] * spin[3 * i + 2]);
    }
}
