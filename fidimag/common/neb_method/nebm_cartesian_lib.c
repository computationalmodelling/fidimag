#include "nebm_cartesian_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_distance_cartesian(double * A, double * B, int n_dofs_image) {

    double A_minus_B[n_dofs_image];
    double distance;

    for(int i = 0; i < n_dofs_image; i++) {
        A_minus_B[i]  = A[i] - B[i];
    }

    distance = compute_norm(A_minus_B, n_dofs_image, n_dofs_image);

    return distance;
}

inline void compute_dYdt(double * m, double * h, double * dm_dt,
                         int * pins,
                         int n_dofs_image
                         ){

    int n_spins = n_dofs_image / 3;
    for(int i = 0; i < n_spins; i++){
       	int j = 3 * i;

        if (pins[i] > 0){
            dm_dt[j] = 0;
		    dm_dt[j + 1] = 0;
		    dm_dt[j + 2] = 0;
		    continue;
		}

        double mm = m[j] * m[j] + m[j + 1] * m[j + 1] + m[j + 2] * m[j + 2];
       	double mh = m[j] * h[j] + m[j + 1] * h[j + 1] + m[j + 2] * h[j + 2];
        //mm.h-mh.m=-mx(mxh)
       	dm_dt[j] = mm * h[j] - mh * m[j];
       	dm_dt[j + 1] = mm * h[j + 1] - mh * m[j + 1];
       	dm_dt[j + 2] = mm * h[j + 2] - mh * m[j + 2];

       	double c = 6 * sqrt(dm_dt[j]     * dm_dt[j]     +
                            dm_dt[j + 1] * dm_dt[j + 1] +
                            dm_dt[j + 2] * dm_dt[j + 2]);

       	dm_dt[j]     += c * (1 - mm) * m[j];
        dm_dt[j + 1] += c * (1 - mm) * m[j + 1];
       	dm_dt[j + 2] += c * (1 - mm) * m[j + 2];
    }

}

void compute_dYdt_C(double * y, double * G, double * dYdt, int * pins, 
                    int n_images, int n_dofs_image) {

	for(int i = 1; i < n_images - 1; i++){

		int j = i * n_dofs_image;

		compute_dYdt(&y[j], &G[j], &dYdt[j], &pins[0], n_dofs_image);

        }

    return;

}
