#include<cmath>
#include "c_vectormath.h"


double cross_x(double a0, double a1, double a2, double b0, double b1, double b2) {
    return a1*b2 - a2*b1;
}

double cross_y(double a0, double a1, double a2, double b0, double b1, double b2) {
    return a2*b0 - a0*b2;
}

double cross_z(double a0, double a1, double a2, double b0, double b1, double b2) {
    return a0*b1 - a1*b0;
}

void normalise(double *a, int n) {
    /* Normalise the *a array, whose length is n (3 * number of nodes in
     * cartesian, and 2 * number of nodes in spherical) To do this we compute
     * the length of *a :
     *
     *      SQRT[ a[0] ** 2 + a[1] ** 2 + ... ]
     *
     *  and divide every *a entry by that length
     */

    double length;

    length = compute_norm(a, n);

    if (length > 0){
        length = 1.0 / length;
    }

    for(int i = 0; i < n; i++){
        a[i] *= length;
    }
}

void normalise(double *m, int *pins, int n) {
    // Normalises n three vectors
    int i, j, k;
    double mdotm;
    for (int id = 0; id < n; id++) {
        i = 3*id;
        j = i + 1;
        k = j + 1;

        if (pins[id]>0) continue;

        mdotm = std::sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
        if(mdotm > 0) {
            mdotm = 1 / mdotm;
            m[i] *= mdotm;
            m[j] *= mdotm;
            m[k] *= mdotm;
        }
    }
}




void cross(double * output, double * A, double * B){
    /* Cross product between arrays A and B, assuming they
     * have length = 3
     * Every resulting component is stored in the *output array
     */

    output[0] = A[1] * B[2] - A[2] * B[1];
    output[1] = A[2] * B[0] - A[0] * B[2];
    output[2] = A[0] * B[1] - A[1] * B[0];
}

double dot(double * A, double * B, int n){
    /* Dot product between arrays A and B , assuming they
     * have length n */

    double dotp = 0;
    for(int i = 0; i < n; i++){
        dotp += A[i] * B[i];
    }

    return dotp;
}

double compute_norm(double * a, int n) {
    /* Compute the norm of an array *a
     *
     * ARGUMENTS:

     * n        :: length of a

     */

    double norm = 0;
    for(int i = 0; i < n; i++) {
        norm += a[i] * a[i];
    }
    norm = std::sqrt(norm);
    return norm;
}
