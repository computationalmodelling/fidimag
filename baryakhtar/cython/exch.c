#include "baryakhtar_clib.h"

void compute_laplace_m(double *m, double *field, double A, double dx, double dy, double dz,
		int nx, int ny, int nz){
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;

	double ax = A/(dx*dx), ay = A/(dy*dy), az = A/(dz*dz);

	#pragma omp parallel for
	for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

            	int id =0;
                int index = nyz * i + nz * j + k;

                double fx=0,fy=0,fz=0;

                if (k > 0) {
                    id = index - 1;
                    fx += az * (m[id] - m[index]);
                    fy += az * (m[id + n1]-m[index+n1]);
                    fz += az * (m[id + n2]-m[index+n2]);
                }

                if (j > 0 ) {
                    id = index - nz;

                    fx += ay * (m[id] - m[index]);
                    fy += ay * (m[id + n1]-m[index+n1]);
                    fz += ay * (m[id + n2]-m[index+n2]);
                }

                if (i > 0) {
                    id = index - nyz;
                    fx += ax * (m[id] - m[index]);
                    fy += ax * (m[id + n1]-m[index+n1]);
                    fz += ax * (m[id + n2]-m[index+n2]);
                }

                if (i < nx - 1) {
                    id = index + nyz;
                    fx += ax * (m[id] - m[index]);
                    fy += ax * (m[id + n1]-m[index+n1]);
                    fz += ax * (m[id + n2]-m[index+n2]);
                }

                if (j < ny - 1 ) {
                    id = index + nz;
                    fx += ay * (m[id] - m[index]);
                    fy += ay * (m[id + n1]-m[index+n1]);
                    fz += ay * (m[id + n2]-m[index+n2]);
                }

                if (k < nz - 1) {
                    id = index + 1;
                    fx += az * (m[id] - m[index]);
                    fy += az * (m[id + n1]-m[index+n1]);
                    fz += az * (m[id + n2]-m[index+n2]);
                }

                field[index] = fx;
                field[index + n1] = fy;
                field[index + n2] = fz;
            }
        }

	}
}


void compute_exch_field_baryakhtar(double *m, double *field, double Me,
						double chi_inv, double A, double dx, double dy, double dz,
						int nx, int ny, int nz) {


	double reduced_A = 2*A/(MU0*Me);

	compute_laplace_m(m, field, reduced_A, dx, dy, dz, nx, ny, nz);

	int n = nx*ny*nz;
	double relax = Me*chi_inv/2.0;

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int j = i+n;
		int k = j+n;

		double mm = m[i]*m[i] + m[j]*m[j] + m[k]*m[k];

		field[i] += relax * (1-mm)* m[i];
		field[j] += relax * (1-mm)* m[j];
		field[k] += relax * (1-mm)* m[k];

	}
}
