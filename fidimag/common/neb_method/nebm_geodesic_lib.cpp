#include "nebm_geodesic_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_geodesic_Vincenty(double * A, double * B,
                                 int n_dofs_image,
                                 int * material,
                                 int n_dofs_image_material
                                 ) {

    /* Compute the Geodesic distance between two images: A and B,  of an energy
     * band. For this, we use Vincenty's formula, which is defined in
     * [Computer Physics Communications 196 (2015) 335–347]
     *
     * INPUTS:
     *
     * A, B         :: Arrays in Cartesian coordinates, i.e.
     *                      A = [ mx_0 my_0 mz_0 mx_1 my_1 mz_1 mx_2...]
     *                 where mx_i is the x component of the i-th spin of image
     *                 A
     *
     * n_dofs_image :: Number of degrees of freedom (total number of spin
     *                 components) per image, which is basically the length of
     *                 the A or B arrays. Since we assume Cartesian
     *                 coordinates, the number of spins is just:
     *                      n_dofs_image / 3
     * 
     * material     :: An array of size (3 * n_spins), containing only 1 and 0s
     *                 For every spin, its value is repeated 3 times (the number
     *                 of dofs). The values of material is 1 where mu_s or M_s
     *                 are larger than zero.
     *                 This array is not necessary here since sites without
     *                 material do not change their magnetisation and thus they
     *                 do not contribute to the Geodesic distance magnitude 
     *
     * Vicenty's formula is defined by computing the distance between
     * corresponding spins of the images A and B. So, if we have m(A)_i as
     * the i-th spin of the A image, then the geodesic distance L is defined
     * as:
     *                _____________________________________________________
     *               /
     *      L = /\  / l(A, B)_0 ^ 2 + l(A, B)_0 ^ 2 + ... + l(A, B)_(N) ^ 2
     *            \/
     *
     *  where N = (n_dofs_image - 1) and
     *
     *                                 ->       ->        ->       ->
     *          l(A, B)_i = arctan2( | m(A)_i X m(B)_i |, m(A)_i o m(B)_i )
     *                                                     (dot product)
     */

    double A_cross_B[3];
    double A_cross_B_norm;
    double A_dot_B;
    int spin_i;
    double geo_dist = 0;
    double distance = 0;
    int n_spins = n_dofs_image / 3;

    // For every spin we will compute the corresponding cross and dot products
    // of A and B. The i-th spin components start at the 3 * i position
    // in the arrays
    // We do not sum sites without
    for(int i = 0; i < n_spins; i++){
        // We only need to know if the spin is in a site with material (Ms, mu_s > 0)
        // so we skip the material for m_y and m_z which are the same as
        // the material for m_x
        if (material[3 * i] > 0) {
            spin_i = 3 * i;

            cross_product(A_cross_B, &A[spin_i], &B[spin_i]);
            A_cross_B_norm = compute_norm(A_cross_B, 3);
            A_dot_B = dot_product(&A[spin_i], &B[spin_i], 3);
            geo_dist = atan2(A_cross_B_norm, A_dot_B);

            distance += geo_dist * geo_dist;
        }
    }

    distance = sqrt(distance);

    return distance;
}

double compute_geodesic_GreatCircle(double * A, double * B,
                                 int n_dofs_image,
                                 int * material,
                                 int n_dofs_image_material
                                 ) {

    /* Compute the Geodesic distance between two images: A and B,  of an energy
     * band. 
     *
     * INPUTS:
     *
     * A, B         :: Arrays in Cartesian coordinates, i.e.
     *                      A = [ mx_0 my_0 mz_0 mx_1 my_1 mz_1 mx_2...]
     *                 where mx_i is the x component of the i-th spin of image
     *                 A
     *
     * n_dofs_image :: Number of degrees of freedom (total number of spin
     *                 components) per image, which is basically the length of
     *                 the A or B arrays. Since we assume Cartesian
     *                 coordinates, the number of spins is just:
     *                      n_dofs_image / 3
     * 
     * material     :: An array of size (3 * n_spins), containing only 1 and 0s
     *                 For every spin, its value is repeated 3 times (the number
     *                 of dofs). The values of material is 1 where mu_s or M_s
     *                 are larger than zero.
     *                 This array is not necessary here since sites without
     *                 material do not change their magnetisation and thus they
     *                 do not contribute to the Geodesic distance magnitude 
     *
     */

    double A_dot_B;
    int spin_i;
    double geo_dist = 0;
    double distance = 0;
    int n_spins = n_dofs_image / 3;

    // For every spin we will compute the corresponding cross and dot products
    // of A and B. The i-th spin components start at the 3 * i position
    // in the arrays
    // We do not sum sites without
    for(int i = 0; i < n_spins; i++){
        // We only need to know if the spin is in a site with material (Ms, mu_s > 0)
        // so we skip the material for m_y and m_z which are the same as
        // the material for m_x
        if (material[3 * i] > 0) {
            spin_i = 3 * i;

            A_dot_B = dot_product(&A[spin_i], &B[spin_i], 3);
            // Prevent nans
            A_dot_B = fmax(-1.0, fmin(1.0, A_dot_B));
            geo_dist = acos(A_dot_B);

            distance += geo_dist * geo_dist;
        }
    }

    distance = sqrt(distance);

    return distance;
}
