#include "baryakhtar_clib.h"

void compute_laplace_m(double *m, double *field, double *Ms, double dx,
		double dy, double dz, int nx, int ny, int nz) {
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;

	double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);

	#pragma omp parallel for
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {

				int id = 0;
				int index = nyz * i + nz * j + k;

				double fx = 0, fy = 0, fz = 0;

				if (Ms[index] == 0.0) {
					field[index] = 0;
					field[index + n1] = 0;
					field[index + n2] = 0;
					continue;
				}

				if (k > 0) {
					id = index - 1;
					if (Ms[id] > 0) {
						fx += az * (m[id] - m[index]);
						fy += az * (m[id + n1] - m[index + n1]);
						fz += az * (m[id + n2] - m[index + n2]);
					}
				}

				if (j > 0) {
					id = index - nz;
					if (Ms[id] > 0) {
						fx += ay * (m[id] - m[index]);
						fy += ay * (m[id + n1] - m[index + n1]);
						fz += ay * (m[id + n2] - m[index + n2]);
					}
				}

				if (i > 0) {
					id = index - nyz;
					if (Ms[id] > 0) {
						fx += ax * (m[id] - m[index]);
						fy += ax * (m[id + n1] - m[index + n1]);
						fz += ax * (m[id + n2] - m[index + n2]);
					}
				}

				if (i < nx - 1) {
					id = index + nyz;
					if (Ms[id] > 0) {
						fx += ax * (m[id] - m[index]);
						fy += ax * (m[id + n1] - m[index + n1]);
						fz += ax * (m[id + n2] - m[index + n2]);
					}
				}

				if (j < ny - 1) {
					id = index + nz;
					if (Ms[id] > 0) {
						fx += ay * (m[id] - m[index]);
						fy += ay * (m[id + n1] - m[index + n1]);
						fz += ay * (m[id + n2] - m[index + n2]);
					}
				}

				if (k < nz - 1) {
					id = index + 1;
					if (Ms[id] > 0) {
						fx += az * (m[id] - m[index]);
						fy += az * (m[id + n1] - m[index + n1]);
						fz += az * (m[id + n2] - m[index + n2]);
					}
				}

				field[index] = fx;
				field[index + n1] = fy;
				field[index + n2] = fz;
			}
		}

	}
}

void compute_relaxation_field_c(double *m, double *field, double *Ms,
		double chi_inv, int n) {

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int j = i + n;
		int k = j + n;

		double relax = Ms[i] * chi_inv / 2.0;
		double mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];

		field[i] = relax * (1 - mm) * m[i];
		field[j] = relax * (1 - mm) * m[j];
		field[k] = relax * (1 - mm) * m[k];

	}
}

void compute_perp_field_c(double *m, double *field, double *field_p, int n) {

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int j = i + n;
		int k = j + n;

		double mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
		double hm = field[i] * m[i] + field[j] * m[j] + field[k] * m[k];

		field_p[i] = field[i] - hm * m[i] / mm;
		field_p[j] = field[j] - hm * m[j] / mm;
		field_p[k] = field[k] - hm * m[k] / mm;

	}
}
