#include "nebm_spherical_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_norm_spherical(double *a, int n, int scale) {
    /* Compute the norm of an array *a. *a is assumed to have spherical
     * coordinates as:
           
     *      a = [ theta_0 phi_0 theta_1 phi_1 ... ]
     
     * This function is necessary to either normalise a vector or compute the
     * distance between two images (see below)
     *
     * To compute the norm, we redefine the angles, so they always lie on the
     * [-PI, PI] range, rather than [0, 2PI] for phi. 
     *
     * ARGUMENTS:
     
     * n        :: length of a
     *
     * scale    :: If scale is zero, we do not scale the norm
     
     * Notice that, when computing the DISTANCE between two images, Y_i, Y_i+1
     * for example, we just do the difference (Y_i - Y_(i+1)) and then compute
     * its norm (using this function) but RESCALING the length by the number
     * of componentsm, thus we use scale=n
     *
     */

    double norm = 0;
    double scale_db = (double) scale;

    for(int i = 0; i < n; i++){

        if (a[i] > WIDE_PI){
            a[i] = 2 * WIDE_PI - a[i];
        } else if(a[i] < -WIDE_PI){
            a[i] += 2 * WIDE_PI;
        }

        norm += a[i] * a[i];
    }

    if (scale == 0) norm = sqrt(norm);
    else norm = sqrt(norm) / scale_db;

    return norm;
}

void normalise_spherical(double *a, int n){

    /* Normalise the *a array, whose length is n (3 * number of nodes in
     * cartesian, and 2 * number of nodes in spherical) To do this we compute
     * the length of *a :
     *
     *      SQRT[ a[0] ** 2 + a[1] ** 2 + ... ]
     *
     *  and divide every *a entry by that length
     */

    double length;

    length = compute_norm_spherical(a, n, 0);

    if (length > 0){
        length = 1.0 / length;
    }

    for(int i = 0; i < n; i++){
        a[i] *= length;
    }
}

/* ------------------------------------------------------------------------- */

double compute_distance_spherical(double * A, double * B, int n) {

    double A_minus_B[n];
    double distance;

    for(int i = 0; i < n; i++) {
        A_minus_B[i]  = A[i] - B[i];
    }

    distance = compute_norm_spherical(A_minus_B, n, n);

    return distance;
}
