#include "micro_clib.h"

void compute_exch_field_micro(double *m, double *field, double *energy,
		double *Ms_inv, double A, double dx, double dy, double dz, int nx,
		int ny, int nz, int xperiodic, int yperiodic) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;

	double ax = 2 * A / (dx * dx), ay = 2 * A / (dy * dy), az = 2 * A / (dz * dz);

	#pragma omp parallel for private(i, j, k, index, id)
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {

				index = nyz * i + nz * j + k;

				if (Ms_inv[index] == 0.0) {
					field[index] = 0;
					field[index + n1] = 0;
					field[index + n2] = 0;
					energy[index] = 0;
					continue;
				}

				double fx = 0, fy = 0, fz = 0;

				if (k > 0) {
					id = index - 1;
					//in general this is unnecessary for normal LLG equation,
					//but for LLbar this is very important.
					if (Ms_inv[id] > 0) {
						fx += az * (m[id] - m[index]);
						fy += az * (m[id + n1] - m[index + n1]);
						fz += az * (m[id + n2] - m[index + n2]);
					}
				}

				if (j > 0 || yperiodic) {
					id = index - nz;
					if (j == 0) {
						id += nyz;
					}
					if (Ms_inv[id] > 0) {
						fx += ay * (m[id] - m[index]);
						fy += ay * (m[id + n1] - m[index + n1]);
						fz += ay * (m[id + n2] - m[index + n2]);
					}
				}

				if (i > 0 || xperiodic) {
					id = index - nyz;
					if (i == 0) {
						id += n1;
					}
					if (Ms_inv[id] > 0) {
						fx += ax * (m[id] - m[index]);
						fy += ax * (m[id + n1] - m[index + n1]);
						fz += ax * (m[id + n2] - m[index + n2]);
					}
				}

				if (i < nx - 1 || xperiodic) {
					id = index + nyz;
					if (i == nx - 1) {
						id -= n1;
					}
					if (Ms_inv[id] > 0) {
						fx += ax * (m[id] - m[index]);
						fy += ax * (m[id + n1] - m[index + n1]);
						fz += ax * (m[id + n2] - m[index + n2]);
					}
				}

				if (j < ny - 1 || yperiodic) {
					id = index + nz;
					if (j == ny - 1) {
						id -= nyz;
					}
					if (Ms_inv[id] > 0) {
						fx += ay * (m[id] - m[index]);
						fy += ay * (m[id + n1] - m[index + n1]);
						fz += ay * (m[id + n2] - m[index + n2]);
					}
				}

				if (k < nz - 1) {
					id = index + 1;
					if (Ms_inv[id] > 0) {
						fx += az * (m[id] - m[index]);
						fy += az * (m[id + n1] - m[index + n1]);
						fz += az * (m[id + n2] - m[index + n2]);
					}
				}

				energy[index] = -0.5 * (fx * m[index] + fy * m[index + n1] + fz * m[index + n2]);

				field[index] = fx * Ms_inv[index] * MU0_INV;
				field[index + n1] = fy * Ms_inv[index] * MU0_INV;
				field[index + n2] = fz * Ms_inv[index] * MU0_INV;

			}
		}
	}
}

