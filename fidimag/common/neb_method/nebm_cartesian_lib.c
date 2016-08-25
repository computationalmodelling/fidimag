#include "nebm_cartesian_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_distance_cartesian(double * A, double * B, int n) {

    double A_minus_B[n];
    double distance;

    for(int i = 0; i < n; i++) {
        A_minus_B[i]  = A[i] - B[i];
    }

    distance = compute_norm(A_minus_B, n, n);

    return distance;
}
