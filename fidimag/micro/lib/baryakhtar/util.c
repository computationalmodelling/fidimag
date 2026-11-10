#include "baryakhtar_clib.h"

/* Vector fields are x0, y0, z0, x1, y1, z1, ... and the cells are numbered with x running fastest, as everywhere else. */
void compute_laplace_m(double *m, double *field, double *Ms, double dx,
		double dy, double dz, int nx, int ny, int nz) {
	const int nxy = nx * ny;
	const double ax = 1.0 / (dx * dx), ay = 1.0 / (dy * dy), az = 1.0 / (dz * dz);

	#pragma omp parallel for
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				const int index = k * nxy + j * nx + i;
				double f[3] = {0, 0, 0};

				if (Ms[index] > 0.0) {
					// (a cell without material is not a neighbour of anything)
					int ids[6], count = 0;
					double a[6];
					if (i > 0) { ids[count] = index - 1; a[count++] = ax; }
					if (i < nx - 1) { ids[count] = index + 1; a[count++] = ax; }
					if (j > 0) { ids[count] = index - nx; a[count++] = ay; }
					if (j < ny - 1) { ids[count] = index + nx; a[count++] = ay; }
					if (k > 0) { ids[count] = index - nxy; a[count++] = az; }
					if (k < nz - 1) { ids[count] = index + nxy; a[count++] = az; }
					for (int q = 0; q < count; q++) {
						if (Ms[ids[q]] > 0.0)
							for (int c = 0; c < 3; c++) f[c] += a[q] * (m[3 * ids[q] + c] - m[3 * index + c]);
					}
				}

				for (int c = 0; c < 3; c++) field[3 * index + c] = f[c];
			}
		}
	}
}

void compute_relaxation_field_c(double *m, double *field, double *Ms,
		double chi_inv, int n) {

	#pragma omp parallel for
	for (int id = 0; id < n; id++) {
		const double *s = m + 3 * id;
		const double relax = Ms[id] * chi_inv / 2.0;
		const double mm = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];

		for (int c = 0; c < 3; c++) field[3 * id + c] = relax * (1 - mm) * s[c];
	}
}

void compute_perp_field_c(double *m, double *field, double *field_p, int n) {

	#pragma omp parallel for
	for (int id = 0; id < n; id++) {
		const double *s = m + 3 * id;
		const double *h = field + 3 * id;

		const double mm = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
		const double hm = h[0] * s[0] + h[1] * s[1] + h[2] * s[2];

		// (a cell with no spin has no perpendicular field)
		for (int c = 0; c < 3; c++) field_p[3 * id + c] = mm > 0.0 ? h[c] - hm * s[c] / mm : 0.0;
	}
}
