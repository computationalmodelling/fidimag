#include "nebm_geodesic_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_distance_geodesic(double * A, double * B, int n_dofs_image){
    /* We will use Vicenty's formula 
     *
     * A, B         :: Arrays in Cartesian coordinates
     *
     */

    int i, j;
    int spin_i;
    double A_cross_B[3];
    double A_cross_B_norm;
    double A_dot_B;
    double distance = 0;
    int n_spins = n_dofs_image / 3;

    for(int i = 0; i < n_spins; i++){
        spin_i = 3 * i;
        cross_product(A_cross_B, &A[spin_i], &B[spin_i]);
        A_cross_B_norm = compute_norm(A_cross_B, 3, 0);

        A_dot_B = dot_product(&A[spin_i], &B[spin_i], 3);

        distance += atan2(A_cross_B_norm, A_dot_B);
    }
    
    distance = sqrt(distance);

    return distance;
}

