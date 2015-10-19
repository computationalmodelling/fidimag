#include "micro_clib.h"

void compute_exch_field_micro(double *m, double *field, double *energy,
			      double *Ms_inv, double A, double dx, double dy, double dz,
                  int nxyz, int *ngbs) {

	double ax = 2 * A / (dx * dx);
    double ay = 2 * A / (dy * dy);
    double az = 2 * A / (dz * dz);

	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {
	    double fx = 0, fy = 0, fz = 0;
	    int idm = 0;
	    int idn = 6 * i; // index for the neighbours

	    if (Ms_inv[i] == 0.0){
	        field[3 * i] = 0;
	        field[3 * i + 1] = 0;
	        field[3 * i + 2] = 0;
	        continue;
	    }

        for (int j = 0; j < 6; j++) {
	        if (ngbs[idn + j] >= 0) {
	            idm = 3 * ngbs[idn + j];
                // if (Ms_inv[i] > 0){

                if (j == 0 || j == 1) {
                    fx += ax * (m[idm]     - m[3 * i]);
                    fy += ax * (m[idm + 1] - m[3 * i + 1]);
                    fz += ax * (m[idm + 2] - m[3 * i + 2]);
                }

                else if (j == 2 || j == 3) {
                    fx += ay * (m[idm]     - m[3  * i]);
                    fy += ay * (m[idm + 1] - m[3 * i + 1]);
                    fz += ay * (m[idm + 2] - m[3 * i + 2]);
                }

                else if (j == 4 || j == 5) {
                    fx += az * (m[idm]     - m[3 * i]);
                    fy += az * (m[idm + 1] - m[3 * i + 1]);
                    fz += az * (m[idm + 2] - m[3 * i + 2]);
                }

                else {
                    printf("Passing");
                    continue; }
                //}
            }
        }

        energy[i] = -0.5 * (fx * m[3 * i] + fy * m[3 * i + 1]
                            + fz * m[3 * i + 2]);

        // printf("energy=%f", energy[i]);

        field[3 * i]     = fx * Ms_inv[i] * MU0_INV;
        field[3 * i + 1] = fy * Ms_inv[i] * MU0_INV;
        field[3 * i + 2] = fz * Ms_inv[i] * MU0_INV;
    }
}
