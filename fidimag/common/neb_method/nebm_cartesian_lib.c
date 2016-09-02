#include "nebm_cartesian_lib.h"
#include "nebm_lib.h"
#include "math.h"

double compute_distance_cartesian(double * A, double * B, int n_dofs_image,
                                  int * material, int n_dofs_image_material
                                  ) {

    /* Compute the distance between two images, A and B, discarding the sites
     * without material
     *
     * We still consider pinned sites because they contribute to the scale
     * factor for the distance
     */

    double A_minus_B[n_dofs_image_material];
    double distance;
    int j = 0;

    for(int i = 0; i < n_dofs_image; i++) {
        if (material[i] > 0) {
            A_minus_B[j]  = A[i] - B[i];
            j += 1;
        }
    }

    distance = compute_norm(A_minus_B, n_dofs_image_material,
                            n_dofs_image_material);

    return distance;
}

inline void compute_dYdt(double * m, double * h, double * dYdt,
                         int * pins,
                         int n_dofs_image
                         ){

    int n_spins = n_dofs_image / 3;
    for(int i = 0; i < n_spins; i++){
       	int j = 3 * i;

        if (pins[i] > 0){
            dYdt[j] = 0;
		    dYdt[j + 1] = 0;
		    dYdt[j + 2] = 0;
		    continue;
		}

        double mm = m[j] * m[j] + m[j + 1] * m[j + 1] + m[j + 2] * m[j + 2];
       	double mh = m[j] * h[j] + m[j + 1] * h[j + 1] + m[j + 2] * h[j + 2];
        //mm.h-mh.m=-mx(mxh)
       	dYdt[j] = mm * h[j] - mh * m[j];
       	dYdt[j + 1] = mm * h[j + 1] - mh * m[j + 1];
       	dYdt[j + 2] = mm * h[j + 2] - mh * m[j + 2];

       	double c = 6 * sqrt(dYdt[j]     * dYdt[j]     +
                            dYdt[j + 1] * dYdt[j + 1] +
                            dYdt[j + 2] * dYdt[j + 2]);

       	dYdt[j]     += c * (1 - mm) * m[j];
        dYdt[j + 1] += c * (1 - mm) * m[j + 1];
       	dYdt[j + 2] += c * (1 - mm) * m[j + 2];
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
